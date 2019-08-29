def writeLabProfile(labs_counts, labs_frequencyPerYear, labs_fractionOfSubjects,labs_units, labs_names,
                    labs_stats, labs_aboveBelowNorm, labs_correlatedLabsCoefficients, labs_abscorrelation,
                    labs_correlatedMedsCoefficients, labs_correlatedProceduresCoefficients, 
                    labs_correlatedDiagnosisCoefficients, labs_correlatedPhenotypesCoefficients, 
                    cohort='All', sex='All', race='All', age_low='All', age_high=None,
                    topN=10, correlationCutoff=0.3):
    """Write out Lab Clinical Profile to JSON File and save locally
    
    Keywords:
    Structures from output of calculateAnyProfile(profileType='labs')
    cohort -- short name for cohort, special characters besides hyphens are prohibited (default 'All')
    sex -- specification of whether this is a 'All', 'Male', or 'Female' sex profile (default 'All')
    race -- specification of whether this is 'All', 'White or Caucasian', 'Black or African American', 'Other' race profile (default 'All')
    age_low -- low age range for this profile (default 'All')
    age_high -- high age range for this profile (default None)
    topN -- integer representing the maximum number of correlations to report in the profile, ranked descending (default 10)
    correlationCutoff -- minimum correlation coefficient value to report for whole profile (default 0.3)
    """   
    import os
    import sys
    import sqlalchemy
    import urllib.parse
    import pandas as pd
    import numpy as np
    import getpass
    from dataclasses import dataclass
    from SciServer import Authentication
    from datetime import datetime
    import json
    from fhir_loader import fhir_loader
    from fhirclient.models import clinicalprofile, fhirreference, identifier, codeableconcept, fhirdate, quantity
    import pymssql
    
    # Initialize  profile
    clinicalProfile = clinicalprofile.ClinicalProfile()
    clinicalProfile.resourceType = 'ClinicalProfile'
    
    if sex == 'M':
        sex = 'Male'
    elif sex =='F':
        sex = 'Female'
    
    # Header info
    if (age_low != 'All'):
        clinicalProfile.id = 'jh-labs-'+cohort+'-'+sex+'-'+race+'-'+str(int(age_low))+'-'+str(int(age_high))
        clinicalProfile.identifier  = [identifier.Identifier({'value': 
                                                              'jh-labs-'+cohort+'-'+sex+'-'+race+'-'+
                                                              str(int(age_low))+'-'+str(int(age_high))})]
        clinicalProfile.cohort = fhirreference.FHIRReference({'reference': 
                                                      'Group/jh-labs-'+cohort+'-'+sex+'-'+race+'-'+str(int(age_low))
                                                              +'-'+str(int(age_high))}) 
    else:
        clinicalProfile.id = 'jh-labs-'+cohort+'-'+sex+'-'+race+'-'+str(age_low)
        clinicalProfile.identifier  = [identifier.Identifier({'value': 
                                                              'jh-labs-'+cohort+'-'+sex+'-'+race+'-'+str(age_low)})]
        clinicalProfile.cohort = fhirreference.FHIRReference({'reference': 
                                                      'Group/jh-labs-'+cohort+'-'+sex+'-'+race+'-'+str(age_low)})
    clinicalProfile.status = 'draft'
    clinicalProfile.population = fhirreference.FHIRReference({'reference': 'Group/jh-labs-'+cohort})
     
    clinicalProfile.date = fhirdate.FHIRDate(str(datetime.now()).replace(' ', 'T'))
    clinicalProfile.reporter = fhirreference.FHIRReference({'reference': 'Organization/JHM',
                           'type': 'Organization',
                           'display': 'Johns Hopkins School of Medicine'})
    ## LABS
    labs = list()
    corrmat = (pd.DataFrame(labs_correlatedLabsCoefficients).unstack(level=[0,1]).corr(min_periods=50)
                        .droplevel(level=0).droplevel(level=0,axis=1))
    lab_names = pd.DataFrame({'lab_name':labs_names}).reset_index()
    lab_counts = pd.DataFrame({'lab_counts':labs_counts}).reset_index().rename({'index':'LAB_LOINC'},axis=1)
    lab_info = lab_names.merge(lab_counts, how='inner', on='LAB_LOINC').set_index('LAB_LOINC')

    for thisLab in lab_info.index:
        
        # Check if STDEV is NaN and skip that lab if so
        if np.isnan(float(labs_stats.loc[thisLab]['std'].median())):
            continue
        
        # Build the profile
        thisCPLab = clinicalprofile.ClinicalProfileLab()
#         try:
        thisCPLab.code = [codeableconcept.CodeableConcept(dict(coding=[dict(system='http://loinc.org', 
                                                                            code=thisLab)],
                                                              text=lab_info.loc[thisLab]['lab_name'][0]))]
        thisCPLab.count = int(lab_info.loc[thisLab]['lab_counts'])
        thisCPLab.frequencyPerYear = round(float(labs_frequencyPerYear.loc[thisLab].mean()),3)
        thisCPLab.fractionOfSubjects = round(float(labs_fractionOfSubjects.loc[thisLab].mean()),3)
        thisCPLab.scalarDistribution = clinicalprofile.ClinicalProfileLabScalarDistribution()
        thisCPLab.scalarDistribution.units = quantity.Quantity(dict(unit=str(labs_units.loc[thisLab][0])))
        thisCPLab.scalarDistribution.min = round(float(labs_stats.loc[thisLab]['min'].min()),3)
        thisCPLab.scalarDistribution.max = round(float(labs_stats.loc[thisLab]['max'].max()),3)
        thisCPLab.scalarDistribution.mean = round(float(labs_stats.loc[thisLab]['mean'].mean()),3)
        thisCPLab.scalarDistribution.median = round(float(labs_stats.loc[thisLab]['median'].median()),3)
        thisCPLab.scalarDistribution.stdDev = round(float(labs_stats.loc[thisLab]['std'].median()),3)
        deciles = list()
        for dec in labs_stats.columns[5:]:
            deciles.append(clinicalprofile.ClinicalProfileLabScalarDistributionDecile(
                                                                dict(nth=int(dec), 
                                                                    value=round(labs_stats.loc[thisLab][dec].mean(),3))))
        thisCPLab.scalarDistribution.decile = deciles

        thisCPLab.scalarDistribution.fractionAboveNormal = round(float(labs_aboveBelowNorm.loc[thisLab].aboveNorm.mean()),3)
        thisCPLab.scalarDistribution.fractionBelowNormal = round(float(labs_aboveBelowNorm.loc[thisLab].belowNorm.mean()),3)

        try:
            yearly_vals = dict()
            for year in corrmat.loc[thisLab].index:
                crosstab = corrmat.loc[(thisLab, year)]
                yearly_vals[year] = (crosstab[crosstab.index.get_level_values(level=1).astype('float') == year]
                                             .droplevel(level=1))

            topNcorrs = pd.DataFrame(yearly_vals).apply(np.mean, axis=1).drop(thisLab).nlargest(topN).round(3)

            entries = list()
            for code, corr in topNcorrs.iteritems():
                if  corr <= correlationCutoff:
                    continue
                otherLoinc = [(dict(coding=[dict(system='http://loinc.org', code=code)],
                                                                  text=str(lab_info.loc[code]['lab_name'][0])))]
                entries.append(dict(labcode=otherLoinc, coefficient=corr))

            if not entries:
                print('No correlated Labs for Lab ', thisLab)
            else:
                thisCPLab.scalarDistribution.correlatedLabs = clinicalprofile.ClinicalProfileLabScalarDistributionCorrelatedLabs(
                                                                dict(topn=topN, 
                                                                     entry=entries))
        except:
            print('No correlated Labs for Lab ', thisLab)

        try:
            topNcorrs = (pd.DataFrame(labs_correlatedMedsCoefficients.loc[thisLab].groupby(['JH_INGREDIENT_RXNORM_CODE'])
                                                                                .Relative_Counts.mean())
                                                                                .Relative_Counts.nlargest(topN).round(3))
            entries = list()
            for code, corr in topNcorrs.iteritems():
                if  corr <= correlationCutoff:
                    continue
                otherRX = [dict(medicationCodeableConcept=dict(coding=
                    [dict(system='http://www.nlm.nih.gov/research/umls/rxnorm/', code=code)]))]
                entries.append(dict(meds=otherRX, coefficient=corr))

            if not entries:
                print('No correlated Meds for Lab ', thisLab)
            else:
                thisCPLab.scalarDistribution.correlatedMedications = clinicalprofile.\
                                        ClinicalProfileLabScalarDistributionCorrelatedMedications(
                                                                        dict(topn=topN, 
                                                                          entry=entries))
        except:
            print('No correlated Meds for Lab ', thisLab)

        try:
            topNcorrs = (pd.DataFrame(labs_correlatedDiagnosisCoefficients.loc[thisLab].groupby(['DX'])
                                                                                .Relative_Counts.mean())
                                                                                .Relative_Counts.nlargest(topN).round(3))
            entries = list()
            for code, corr in topNcorrs.iteritems():
                if  corr <= correlationCutoff:
                    continue
                otherDX = (dict(coding=[dict(system='http://www.icd10data.com/', code=code)]))
                entries.append(dict(code=otherDX, coefficient=corr))

            if not entries:
                print('No correlated Diagnoses for Lab ', thisLab)
            else:
                thisCPLab.scalarDistribution.correlatedDiagnoses = clinicalprofile.\
                                                            ClinicalProfileLabScalarDistributionCorrelatedDiagnoses(
                                                                    dict(topn=topN, 
                                                                      entry=entries))
        except:
            print('No correlated Diagnoses for Lab ', thisLab)

        try:      
            topNcorrs = (pd.DataFrame(labs_correlatedProceduresCoefficients.loc[thisLab].groupby(['RAW_PX'])
                                                                                .Relative_Counts.mean())
                                                                                .Relative_Counts.nlargest(topN).round(3))
            entries = list()
            for code, corr in topNcorrs.iteritems():
                if  corr <= correlationCutoff:
                    continue
                otherProc = [(dict(coding=[dict(system='http://www.ama-assn.org/practice-management/cpt', code=code)]))]
                entries.append(dict(code=otherProc, coefficient=corr))

            if not entries:
                print('No correlated Procedures for Lab ', thisLab)
            else:
                thisCPLab.scalarDistribution.correlatedProcedures = clinicalprofile.\
                                                            ClinicalProfileLabScalarDistributionCorrelatedProcedures(
                                                                    dict(topn=topN, 
                                                                      entry=entries))
        except:
            print('No correlated Procedures for Lab ', thisLab)

        try:      
            topNcorrs = (pd.DataFrame(labs_correlatedPhenotypesCoefficients.loc[thisLab].groupby(['HPO'])
                                                                                .Relative_Counts.mean())
                                                                                .Relative_Counts.nlargest(topN).round(3))
            entries = list()
            for code, corr in topNcorrs.iteritems():
                if  corr <= correlationCutoff:
                    continue
                otherHPO = (dict(coding=[dict(system='http://hpo.jax.org/app/', code=code)]))
                entries.append(dict(code=otherHPO, coefficient=corr))

            if not entries:
                print('No correlated Phenotypes for Lab ', thisLab)
            else:
                thisCPLab.scalarDistribution.correlatedPhenotypes = clinicalprofile.\
                                                            ClinicalProfileLabScalarDistributionCorrelatedPhenotypes(
                                                                    dict(topn=topN, 
                                                                      entry=entries))
        except:
            print('No correlated Phenotypes for Lab ', thisLab)

        labs.append(thisCPLab)
        
#         except:
#             print('This lab did not work ', thisLab)
        
    clinicalProfile.lab = labs

    if age_high != None:
        filename = cohort+'_resources/jh-labs-'+cohort+'-'+sex+'-'+race+'-'+str(int(age_low))+'-'+str(int(age_high))+'.json'
    else:
        filename = cohort+'_resources/jh-labs-'+cohort+'-'+sex+'-'+race+'-'+str(age_low)+'.json'
        
    with open(filename, 'w') as outfile:
        json.dump(clinicalProfile.as_json(), outfile, indent=4)
    
    del(clinicalProfile)
    return print('Write to '+ filename + ' successful')