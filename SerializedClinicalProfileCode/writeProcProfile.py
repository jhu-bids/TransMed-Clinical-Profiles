def writeProcProfile(procedures_code, procedures_count, procedures_frequencyPerYear, procedures_fractionOfSubjects,
                    procs_correlatedLabsCoefficients, procs_correlatedDiagsCoefficients, procs_correlatedMedsCoefficients,
                    procs_correlatedProceduresCoefficients, procs_correlatedPhenotypesCoefficients,
                     cohort='All', sex='All', race='All', age_low='All', age_high=None,
                    topN=10, correlationCutoff=0.3):
    """Write out Procedures Clinical Profile to JSON File and save locally
    
    Keywords:
    Structures from output of calculateAnyProfile(profileType='procedures')
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
        clinicalProfile.id = 'jh-procedures-'+cohort+'-'+sex+'-'+race+'-'+str(int(age_low))+'-'+str(int(age_high))
        clinicalProfile.identifier  = [identifier.Identifier({'value': 
                                                              'jh-procedures-'+cohort+'-'+sex+'-'+race+'-'+str(int(age_low))+'-'+str(int(age_high))})]
        clinicalProfile.cohort = fhirreference.FHIRReference({'reference': 
                                                      'Group/jh-procedures-'+cohort+'-'+sex+'-'+race+'-'+str(int(age_low))+'-'+str(int(age_high))}) 
    else:
        clinicalProfile.id = 'jh-procedures-'+cohort+'-'+sex+'-'+race+'-'+str(age_low)
        clinicalProfile.identifier  = [identifier.Identifier({'value': 
                                                              'jh-procedures-'+cohort+'-'+sex+'-'+race+'-'+str(age_low)})]
        clinicalProfile.cohort = fhirreference.FHIRReference({'reference': 
                                                      'Group/jh-procedures-'+cohort+'-'+sex+'-'+race+'-'+str(age_low)})
    clinicalProfile.status = 'draft'
    clinicalProfile.population = fhirreference.FHIRReference({'reference': 'Group/jh-procedures-'+cohort})
     
    clinicalProfile.date = fhirdate.FHIRDate(str(datetime.now()).replace(' ', 'T'))
    clinicalProfile.reporter = fhirreference.FHIRReference({'reference': 'Organization/JHM',
                           'type': 'Organization',
                           'display': 'Johns Hopkins School of Medicine'})
    
    procs = list()
    for thisProc in procedures_code:
        thisCPProc = clinicalprofile.ClinicalProfileProcedure()
        try:
            thisCPProc.code = [codeableconcept.CodeableConcept(dict(coding=[dict(
                                                            system='http://www.ama-assn.org/practice-management/cpt', 
                                                            code=str(thisProc))]))]

            thisCPProc.frequencyPerYear = round(float(procedures_frequencyPerYear.loc[thisProc].mean()),3)
            thisCPProc.fractionOfSubjects = round(float(procedures_fractionOfSubjects.loc[thisProc].mean()),3)
            
            try:
                topNcorrs = (pd.DataFrame(procs_correlatedLabsCoefficients.loc[thisProc].groupby(['LAB_LOINC'])
                                                                                    .Relative_Counts.mean())
                                                                                    .Relative_Counts.nlargest(topN).round(3))
                entries = list()
                for code, corr in topNcorrs.iteritems():
                    if  corr <= correlationCutoff:
                        continue
                    otherLab = [(dict(coding=[dict(system='http://loinc.org', code=code)]))]
                    entries.append(dict(labcode=otherLab, coefficient=corr))

                if not entries:
                    print('No correlated Labs for Procedure ', thisProc)
                else:
                    thisCPProc.correlatedLabs = clinicalprofile.\
                                    ClinicalProfileLabScalarDistributionCorrelatedLabs(dict(topn=topN, entry=entries))
            except:
                print('No correlated Labs for Procedure ', thisProc)

            try:
                topNcorrs = (pd.DataFrame(procs_correlatedDiagsCoefficients.loc[thisProc].groupby(['DX'])
                                                                                    .Relative_Counts.mean())
                                                                                    .Relative_Counts.nlargest(topN).round(3))
                entries = list()
                for code, corr in topNcorrs.iteritems():
                    if  corr <= correlationCutoff:
                        continue
                    otherDX = (dict(coding=[dict(system='http://www.icd10data.com/', code=code)]))
                    entries.append(dict(code=otherDX, coefficient=corr))

                if not entries:
                    print('No correlated Diagnoses for Procedure ', thisProc)
                else:
                    thisCPProc.correlatedDiagnoses = clinicalprofile.ClinicalProfileLabScalarDistributionCorrelatedDiagnoses(dict
                                                                                                        (topn=topN, 
                                                                                                         entry=entries))
            except Exception as e:
                print(e)
                print('No correlated DX for Procedure ', thisProc)

            try:
                topNcorrs = (pd.DataFrame(procs_correlatedProceduresCoefficients.loc[thisProc].groupby(['RAW_PX'])
                                                                                    .Relative_Counts.mean())
                                                                                    .Relative_Counts.nlargest(topN).round(3))
                entries = list()
                for code, corr in topNcorrs.iteritems():
                    if  corr <= correlationCutoff:
                        continue
                    otherProc = [(dict(coding=[dict(system='http://www.ama-assn.org/practice-management/cpt', code=code)]))]
                    entries.append(dict(code=otherProc, coefficient=corr))

                if not entries:
                    print('No correlated Procedures for Procedure ', thisProc)
                else:
                    thisCPProc.correlatedProcedures = clinicalprofile.ClinicalProfileLabScalarDistributionCorrelatedProcedures(dict
                                                                                                        (topn=topN, 
                                                                                                         entry=entries))
            except Exception as e:
                print(e)
                print('No correlated Procedures for Procedure ', thisProc)

            try:
                topNcorrs = (pd.DataFrame(procs_correlatedMedsCoefficients.loc[thisProc].groupby(['JH_INGREDIENT_RXNORM_CODE'])
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
                    print('No correlated Meds for Procedure ', thisProc)
                else:
                    thisCPProc.correlatedMedications = clinicalprofile.ClinicalProfileLabScalarDistributionCorrelatedMedications(dict
                                                                                                        (topn=topN, 
                                                                                                         entry=entries))
            except Exception as e:
                print(e)
                print('No correlated Meds for Procedure ', thisProc)

            try:      
                topNcorrs = (pd.DataFrame(procs_correlatedPhenotypesCoefficients.loc[thisProc].groupby(['HPO'])
                                                                                    .Relative_Counts.mean())
                                                                                    .Relative_Counts.nlargest(topN).round(3))
                entries = list()
                for code, corr in topNcorrs.iteritems():
                    if  corr <= correlationCutoff:
                        continue
                    otherHPO = (dict(coding=[dict(system='http://hpo.jax.org/app/', code=code)]))
                    entries.append(dict(code=otherHPO, coefficient=corr))

                if not entries:
                    print('No correlated Phenotypes for Procedure ', thisProc)
                else:
                    thisCPProc.correlatedPhenotypes = clinicalprofile.ClinicalProfileLabScalarDistributionCorrelatedPhenotypes(
                                                                        dict(topn=topN, 
                                                                          entry=entries))
            except:
                print('No correlated Phenotypes for Procedure ', thisProc)
                
            procs.append(thisCPProc)
            
        except:
            print('This procedure did not work ', thisProc)
        
    clinicalProfile.procedure = procs
    
    if age_high != None:
        filename = cohort+'_resources/jh-procedures-'+cohort+'-'+sex+'-'+race+'-'+str(int(age_low))+'-'+str(int(age_high))+'.json'
    else:
        filename = cohort+'_resources/jh-procedures-'+cohort+'-'+sex+'-'+race+'-'+str(age_low)+'.json'
        
    with open(filename, 'w') as outfile:
        json.dump(clinicalProfile.as_json(), outfile, indent=4)
    
    del(clinicalProfile)
    return print('Write to '+ filename + ' successful')