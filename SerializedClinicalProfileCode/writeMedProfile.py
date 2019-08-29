def writeMedProfile(meds_medication, meds_frequencyPerYear, meds_fractionOfSubjects, 
                    meds_correlatedLabsCoefficients, meds_correlatedMedsCoefficients, 
                    meds_correlatedProceduresCoefficients, meds_correlatedDiagnosisCoefficients,
                    meds_correlatedPhenotypesCoefficients,
                    cohort='All', sex='All', race='All', age_low='All', age_high=None,
                    topN=10, correlationCutoff=0.3):
    """Write out Medication Clinical Profile to JSON File and save locally
    
    Keywords:
    Structures from output of calculateAnyProfile(profileType='medications')
    cohort -- short name for cohort, special characters besides hyphens are prohibited (default 'All')
    sex -- specification of whether this is a 'All', 'Male', or 'Female' sex profile (default 'All')
    race -- specification of whether this is 'All', 'White or Caucasian', 'Black or African American', 'Other'
    race profile (default 'All')
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
    from fhirclient.models import clinicalprofile, fhirreference, identifier, codeableconcept, fhirdate, quantity
    from datetime import datetime
    import json
    from fhir_loader import fhir_loader
    import pymssql
    
    # Initialize  profile
    clinicalProfile = clinicalprofile.ClinicalProfile()
    
    if sex == 'M':
        sex = 'Male'
    elif sex =='F':
        sex = 'Female'
    
    # Header info
    if (age_low != 'All'):
        clinicalProfile.id = ('jh-medications-'+cohort+'-'+sex+'-'+race+'-'
                              +str(int(age_low))+'-'+str(int(age_high)))
        clinicalProfile.identifier  = [identifier.Identifier({'value': 
                                                              ('jh-medications-'+cohort+'-'
                                                              +sex+'-'+race+'-'+str(int(age_low))
                                                              +'-'+str(int(age_high)))})]
        clinicalProfile.cohort = fhirreference.FHIRReference({'reference': 
                                                      ('Group/jh-medications-'+cohort+'-'+sex
                                                              +'-'+race+'-'+str(int(age_low))
                                                              +'-'+str(int(age_high)))}) 
    else:
        clinicalProfile.id = 'jh-medications-'+cohort+'-'+sex+'-'+race+'-'+str(age_low)
        clinicalProfile.identifier  = [identifier.Identifier({'value': 
                                                              'jh-medications-'+cohort+'-'
                                                              +sex+'-'+race+'-'+str(age_low)})]
        clinicalProfile.cohort = fhirreference.FHIRReference({'reference': 
                                                      ('Group/jh-medications-'+cohort+'-'+sex+
                                                              '-'+race+'-'+str(age_low))})
    clinicalProfile.status = 'draft'
    clinicalProfile.population = fhirreference.FHIRReference({'reference': 'Group/jh-medications-'+cohort})
     
    clinicalProfile.date = fhirdate.FHIRDate(str(datetime.now()).replace(' ', 'T'))
    clinicalProfile.reporter = fhirreference.FHIRReference({'reference': 'Organization/JHM',
                           'type': 'Organization',
                           'display': 'Johns Hopkins School of Medicine'})

    meds = list()
    meds_medication = [x for x in meds_medication if str(x) != 'nan']
    for thisMed in meds_medication:
        thisCPMed = clinicalprofile.ClinicalProfileMedication()
        try:
            thisCPMed.medicationCodeableConcept=codeableconcept.CodeableConcept(dict(coding=
                     [dict(system='http://www.nlm.nih.gov/research/umls/rxnorm/', code=str(thisMed))]))
            
            # meds_dosageInfo
#             thisCPMed.dose = quantity.Quantity(dict(unit=str(meds_dosageInfo.loc[thisMed])))
            
            thisCPMed.frequencyPerYear = round(float(meds_frequencyPerYear.loc[thisMed].mean()),3)
            thisCPMed.fractionOfSubjects = round(float(meds_fractionOfSubjects.loc[thisMed].mean()),3)

            try:
                topNcorrs = (pd.DataFrame(meds_correlatedLabsCoefficients.loc[thisMed].groupby(['LAB_LOINC'])
                                                                                    .Relative_Counts.mean())
                                                                                    .Relative_Counts.nlargest(topN).round(3))
                entries = list()
                for code, corr in topNcorrs.iteritems():
                    if  corr <= correlationCutoff:
                        continue
                    otherLab = [(dict(coding=[dict(system='http://loinc.org', code=code)]))]
                    entries.append(dict(labcode=otherLab, coefficient=corr))

                if not entries:
                    print('No correlated Labs for Med ', thisMed)
                else:
                    thisCPMed.correlatedLabs = clinicalprofile.\
                                        ClinicalProfileLabScalarDistributionCorrelatedLabs(dict(topn=topN, entry=entries))
            except:
                print('No correlated Labs for Med ', thisMed)
                
            try:
                topNcorrs = (pd.DataFrame(meds_correlatedDiagnosisCoefficients.loc[thisMed].groupby(['DX'])
                                                                                    .Relative_Counts.mean())
                                                                                    .Relative_Counts.nlargest(topN).round(3))
                entries = list()
                for code, corr in topNcorrs.iteritems():
                    if  corr <= correlationCutoff:
                        continue
                    otherDX = (dict(coding=[dict(system='http://www.icd10data.com/', code=code)]))
                    entries.append(dict(code=otherDX, coefficient=corr))

                if not entries:
                    print('No correlated Diagnoses for Med ', thisMed)
                else:
                    thisCPMed.correlatedDiagnoses = clinicalprofile.ClinicalProfileLabScalarDistributionCorrelatedDiagnoses(dict
                                                                                                        (topn=topN, 
                                                                                                         entry=entries))
            except:
                print('No correlated DX for Med ', thisMed)

            try:
                topNcorrs = (pd.DataFrame(meds_correlatedProceduresCoefficients.loc[thisMed].groupby(['RAW_PX'])
                                                                                    .Relative_Counts.mean())
                                                                                    .Relative_Counts.nlargest(topN).round(3))
                entries = list()
                for code, corr in topNcorrs.iteritems():
                    if  corr <= correlationCutoff:
                        continue
                    otherProc = [(dict(coding=[dict(system='http://www.ama-assn.org/practice-management/cpt', code=code)]))]
                    entries.append(dict(code=otherProc, coefficient=corr))

                if not entries:
                    print('No correlated Procedures for Med ', thisMed)
                else:
                    thisCPMed.correlatedProcedures = clinicalprofile.ClinicalProfileLabScalarDistributionCorrelatedProcedures(dict
                                                                                                        (topn=topN, 
                                                                                                         entry=entries))
            except:
                print('No correlated Procedures for Med ', thisMed)

            try:
                topNcorrs = (pd.DataFrame(meds_correlatedMedsCoefficients.loc[thisMed].groupby(['JH_INGREDIENT_RXNORM_CODE'])
                                                                                    .Relative_Counts.mean())
                                                                                    .Relative_Counts.nlargest(topN).round(3))
                entries = list()
                for code, corr in topNcorrs.iteritems():
                    if  corr <= correlationCutoff:
                        continue
                    otherMed = [dict(medicationCodeableConcept=dict(coding=
                     [dict(system='http://www.nlm.nih.gov/research/umls/rxnorm/', code=code)]))]
                    entries.append(dict(meds=otherMed, coefficient=corr))

                if not entries:
                    print('No correlated Meds for Med ', thisMed)
                else:
                    thisCPMed.correlatedMedications = clinicalprofile.ClinicalProfileLabScalarDistributionCorrelatedMedications(dict
                                                                                                        (topn=topN, 
                                                                                                         entry=entries))
            except:
                print('No correlated Meds for Med ', thisMed)
                
            try:      
                topNcorrs = (pd.DataFrame(meds_correlatedPhenotypesCoefficients.loc[thisMed].groupby(['HPO'])
                                                                                    .Relative_Counts.mean())
                                                                                    .Relative_Counts.nlargest(topN).round(3))
                entries = list()
                for code, corr in topNcorrs.iteritems():
                    if  corr <= correlationCutoff:
                        continue
                    otherHPO = (dict(coding=[dict(system='http://hpo.jax.org/app/', code=code)]))
                    entries.append(dict(code=otherHPO, coefficient=corr))

                if not entries:
                    print('No correlated Phenotypes for Med ', thisMed)
                else:
                    thisCPMed.correlatedPhenotypes = clinicalprofile.ClinicalProfileLabScalarDistributionCorrelatedPhenotypes(
                                                                        dict(topn=topN, 
                                                                          entry=entries))
            except:
                print('No correlated Phenotypes for Med ', thisMed)

            meds.append(thisCPMed)
            
        except Exception as e:
            print(e)
            print('This med did not work ', thisMed)
        
    clinicalProfile.medication = meds
    
    if age_high != None:
        filename = cohort+'_resources/jh-medications-'+cohort+'-'+sex+'-'+race+'-'+str(int(age_low))+'-'+str(int(age_high))+'.json'
    else:
        filename = cohort+'_resources/jh-medications-'+cohort+'-'+sex+'-'+race+'-'+str(age_low)+'.json'
        
    with open(filename, 'w') as outfile:
        json.dump(clinicalProfile.as_json(), outfile, indent=4)
    
    del(clinicalProfile)
    return print('Write to '+ filename + ' successful')