def calculateAnyProfile(profileType, df_labs, df_meds, df_procedures, df_diagnoses, df_phenotypes):
    """Calculate a single profile based on the type provided and data cleaned from getSubdemographicsTables
    
    Arguments:
    profileType -- which individual profile type you would like generated, this will be the category with the header information
    (Options: 'labs', 'medications', 'procedures', 'diagnoses', 'phenotypes')
    
    Keywords:
    df_labs -- labs dataframe returned from getSubdemographicsTables
    df_medications -- medications dataframe returned from getSubdemographicsTables
    df_procedures -- procedures dataframe returned from getSubdemographicsTables
    df_diagnoses -- diagnoses dataframe returned from getSubdemographicsTables
    df_phenotypes -- phenotypes dataframe returned from getSubdemographicsTables
    
    Returns Pythonic structures needed to generate profile in JSON format using the corresponding write profile function
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
    import pymssql
    
    try:
        # Make Labs Profile
        if profileType == 'labs':
            # High Level Info, Scalar Distribution
            labs_counts = df_labs.LAB_LOINC.value_counts()

            grouped_labs = df_labs.groupby(['LAB_LOINC', 'resultYear'])

            labs_frequencyPerYear = (df_labs.groupby(['LAB_LOINC','PATID','resultYear']).PATID.size()
                                            .groupby(['LAB_LOINC','resultYear']).aggregate(np.mean))
            labs_fractionOfSubjects = (np.divide(df_labs.groupby(['LAB_LOINC']).PATID.nunique(),
                                                      df_labs.PATID.nunique()))
            labs_units = df_labs.groupby(['LAB_LOINC']).LOINC_UNIT.unique()
            labs_names = df_labs.groupby(['LAB_LOINC']).LOINC_SHORTNAME.unique()

            def percentile(n):
                def percentile_(x):
                    return x.quantile(n*0.01)
                percentile_.__name__ = '%s' % n
                return percentile_

            labs_stats = (grouped_labs
                       .RESULT_NUM.agg(['min','max', 'mean','median','std',
                                           percentile(10), percentile(20), percentile(30),
                                           percentile(40), percentile(50), percentile(60),
                                           percentile(70), percentile(80), percentile(90)]))

            def fracsAboveBelowNormal(x):
                try:
                    aboveNorm = np.divide(np.sum(x.RESULT_NUM > x.range_high), x.RESULT_NUM.size)
                    belowNorm = np.divide(np.sum(x.RESULT_NUM < x.range_low), x.RESULT_NUM.size)
                    return pd.Series({'aboveNorm':aboveNorm, 'belowNorm':belowNorm})
                except:
                    return pd.Series({'aboveNorm':np.nan, 'belowNorm':np.nan})

            labs_aboveBelowNorm = (grouped_labs.apply(fracsAboveBelowNormal))

            labs_correlatedLabsCoefficients = (df_labs.groupby(['LAB_LOINC','resultYear','PATID'])
                                               .RESULT_NUM.mean())

            labs_abscorrelation = 0

            ## LABS TO MEDICATIONS
            def patientsAboveBelowNormalLabsMeds(x):
                # Get patients above and below normal
                patientsAboveNorm = x.PATID[x.RESULT_NUM > x.range_high].tolist()
                patientsBelowNorm = x.PATID[x.RESULT_NUM < x.range_low].tolist()

                # Get unique patient IDs for above & below normal
                patientsAboveBelowNorm = list(set(patientsAboveNorm + patientsBelowNorm))

                # Link to meds table
                abnormalPatientsMeds = df_meds[df_meds.PATID.isin(patientsAboveBelowNorm) &
                                             (df_meds.startYear == pd.to_datetime(x.RESULT_DATE).dt.year.unique()[0])]

                return pd.Series({'medsAboveBelowNorm': abnormalPatientsMeds.JH_INGREDIENT_RXNORM_CODE.value_counts().index,
                                 'counts': abnormalPatientsMeds.JH_INGREDIENT_RXNORM_CODE.value_counts().values})

            # Need to grab the indices of those with abnormal lab, grab their medications, count and rank them 
            labs_correlatedMedsCoefficients = (grouped_labs.apply(patientsAboveBelowNormalLabsMeds))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for lab in labs_correlatedMedsCoefficients.index:
                thisLabYear = labs_correlatedMedsCoefficients.loc[lab]
                thisLab = lab[0]
                thisYear = lab[1]
                totalCrossTab = np.sum(thisLabYear.counts)
                for medInd in range(len(labs_correlatedMedsCoefficients.loc[lab].medsAboveBelowNorm.values)):
                    mytups.append((thisLabYear.medsAboveBelowNorm.values[medInd], thisLabYear.counts[medInd]/totalCrossTab))
                    multiIndex.append((thisLab, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)
            labs_correlatedMedsCoefficients = (pd.DataFrame.from_records(mytups, columns=['JH_INGREDIENT_RXNORM_CODE','Relative_Counts'],
                                                           index=index))

            ## LABS TO PROCEDURES
            def patientsAboveBelowNormalLabsProcs(x):
                # Get patients above and below normal
                patientsAboveNorm = x.PATID[x.RESULT_NUM > x.range_high].tolist()
                patientsBelowNorm = x.PATID[x.RESULT_NUM < x.range_low].tolist()

                # Get unique patient IDs for above & below normal
                patientsAboveBelowNorm = list(set(patientsAboveNorm + patientsBelowNorm))

                # Link to procs table
                abnormalPatientsProcs = df_procedures[df_procedures.PATID.isin(patientsAboveBelowNorm) &
                                             (df_procedures.encounterYear == pd.to_datetime(x.RESULT_DATE).dt.year.unique()[0])]

                return pd.Series({'procsAboveBelowNorm': abnormalPatientsProcs.RAW_PX.value_counts().index,
                                 'counts': abnormalPatientsProcs.RAW_PX.value_counts().values})

            # Need to grab the indices of those with abnormal lab, grab their medications, count and rank them 
            labs_correlatedProceduresCoefficients = (grouped_labs.apply(patientsAboveBelowNormalLabsProcs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for lab in labs_correlatedProceduresCoefficients.index:
                thisLabYear = labs_correlatedProceduresCoefficients.loc[lab]
                thisLab = lab[0]
                thisYear = lab[1]
                totalCrossTab = np.sum(thisLabYear.counts)
                for procInd in range(len(labs_correlatedProceduresCoefficients.loc[lab].procsAboveBelowNorm.values)):
                    mytups.append((thisLabYear.procsAboveBelowNorm.values[procInd], thisLabYear.counts[procInd]/totalCrossTab))
                    multiIndex.append((thisLab, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)    
            labs_correlatedProceduresCoefficients = (pd.DataFrame.from_records(mytups, columns=['RAW_PX','Relative_Counts'],
                                                                              index=index))

            ## LABS TO DIAGNOSES
            def patientsAboveBelowNormalLabsDiags(x):
                # Get patients above and below normal
                patientsAboveNorm = x.PATID[x.RESULT_NUM > x.range_high].tolist()
                patientsBelowNorm = x.PATID[x.RESULT_NUM < x.range_low].tolist()

                # Get unique patient IDs for above & below normal
                patientsAboveBelowNorm = list(set(patientsAboveNorm + patientsBelowNorm))

                # Link to procs table
                abnormalPatientsDiags = df_diagnoses[df_diagnoses.PATID.isin(patientsAboveBelowNorm) &
                                             (df_diagnoses.admitYear == pd.to_datetime(x.RESULT_DATE).dt.year.unique()[0])]

                return pd.Series({'diagsAboveBelowNorm': abnormalPatientsDiags.DX.value_counts().index,
                                 'counts': abnormalPatientsDiags.DX.value_counts().values})

            # Need to grab the indices of those with abnormal lab, grab their medications, count and rank them 
            labs_correlatedDiagnosisCoefficients = (grouped_labs.apply(patientsAboveBelowNormalLabsDiags))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for lab in labs_correlatedDiagnosisCoefficients.index:
                thisLabYear = labs_correlatedDiagnosisCoefficients.loc[lab]
                thisLab = lab[0]
                thisYear = lab[1]
                totalCrossTab = np.sum(thisLabYear.counts)
                for diagInd in range(len(labs_correlatedDiagnosisCoefficients.loc[lab].diagsAboveBelowNorm.values)):
                    mytups.append((thisLabYear.diagsAboveBelowNorm.values[diagInd], thisLabYear.counts[diagInd]/totalCrossTab))
                    multiIndex.append((thisLab, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            labs_correlatedDiagnosisCoefficients = (pd.DataFrame.from_records(mytups, columns=['DX','Relative_Counts'],
                                                                             index=index))

            ## LABS TO PHENOTYPES
            def patientsAboveBelowNormalLabsHPOs(x):
                # Get patients above and below normal
                patientsAboveNorm = x.PATID[x.RESULT_NUM > x.range_high].tolist()
                patientsBelowNorm = x.PATID[x.RESULT_NUM < x.range_low].tolist()

                # Get unique patient IDs for above & below normal
                patientsAboveBelowNorm = list(set(patientsAboveNorm + patientsBelowNorm))

                # Link to procs table
                abnormalPatientsHPOs = df_phenotypes[df_phenotypes.PATID.isin(patientsAboveBelowNorm) &
                                             (df_phenotypes.admitYear == pd.to_datetime(x.RESULT_DATE).dt.year.unique()[0])]

                return pd.Series({'hposAboveBelowNorm': abnormalPatientsHPOs.HPO.value_counts().index,
                                 'counts': abnormalPatientsHPOs.HPO.value_counts().values})

            # Need to grab the indices of those with abnormal lab, grab their medications, count and rank them 
            labs_correlatedPhenotypesCoefficients = (grouped_labs.apply(patientsAboveBelowNormalLabsHPOs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for lab in labs_correlatedPhenotypesCoefficients.index:
                thisLabYear = labs_correlatedPhenotypesCoefficients.loc[lab]
                thisLab = lab[0]
                thisYear = lab[1]
                totalCrossTab = np.sum(thisLabYear.counts)
                for hpoInd in range(len(labs_correlatedPhenotypesCoefficients.loc[lab].hposAboveBelowNorm.values)):
                    mytups.append((thisLabYear.hposAboveBelowNorm.values[hpoInd], thisLabYear.counts[hpoInd]/totalCrossTab))
                    multiIndex.append((thisLab, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            labs_correlatedPhenotypesCoefficients = (pd.DataFrame.from_records(mytups, columns=['HPO','Relative_Counts'],
                                                                             index=index))

            return (labs_counts, labs_frequencyPerYear, labs_fractionOfSubjects, labs_units, labs_names,
                   labs_stats, labs_aboveBelowNorm, labs_correlatedLabsCoefficients, labs_abscorrelation,
                    labs_correlatedMedsCoefficients, labs_correlatedProceduresCoefficients, labs_correlatedDiagnosisCoefficients,
                   labs_correlatedPhenotypesCoefficients)
            
        # Make Medication Profile
        elif profileType == 'medications':
            meds_medication = df_meds.JH_INGREDIENT_RXNORM_CODE.unique()

            meds_dosageInfo = df_meds.groupby('JH_INGREDIENT_RXNORM_CODE').RX_DOSE_ORDERED.mean()

            meds_frequencyPerYear = (df_meds.groupby(['JH_INGREDIENT_RXNORM_CODE','startYear','PATID']).PATID
                                .count().groupby(['JH_INGREDIENT_RXNORM_CODE','startYear']).mean())

            meds_fractionOfSubjects = (np.divide(df_meds.groupby(['JH_INGREDIENT_RXNORM_CODE']).PATID.nunique(),
                                            df_meds.PATID.nunique()))

            grouped_meds = df_meds.groupby(['JH_INGREDIENT_RXNORM_CODE', 'startYear'])

            #meds_correlatedLabsCoefficients
            def patientsAboveBelowNormalMedsLabs(x):

                patientsWithThisRX = list(set(x.PATID.tolist()))

                # Link to labs table
                abnormalPatientsLabs = df_labs[(df_labs.PATID.isin(patientsWithThisRX)) & 
                                               ((df_labs.RESULT_NUM > df_labs.range_high) | 
                                                (df_labs.RESULT_NUM < df_labs.range_low)) &
                                              (df_labs.resultYear == pd.to_datetime(x.RX_START_DATE).dt.year.unique()[0])]

                return pd.Series({'labsAboveBelowNorm': abnormalPatientsLabs.LAB_LOINC.value_counts().index,
                                 'counts': abnormalPatientsLabs.LAB_LOINC.value_counts().values})

            meds_correlatedLabsCoefficients = (grouped_meds.apply(patientsAboveBelowNormalMedsLabs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for med in meds_correlatedLabsCoefficients.index:
                thisMedYear = meds_correlatedLabsCoefficients.loc[med]
                thisMed = med[0]
                thisYear = med[1]
                totalCrossTab = np.sum(thisMedYear.counts)
                for labInd in range(len(meds_correlatedLabsCoefficients.loc[med].labsAboveBelowNorm.values)):
                    mytups.append((thisMedYear.labsAboveBelowNorm.values[labInd], thisMedYear.counts[labInd]/totalCrossTab))
                    multiIndex.append((thisMed, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            meds_correlatedLabsCoefficients = (pd.DataFrame.from_records(mytups, columns=['LAB_LOINC','Relative_Counts'],
                                                                             index=index))

            #meds_correlatedDiagsCoefficients
            def patientsCrossFreqMedsDiags(x):

                patientsWithThisRX = list(set(x.PATID.tolist()))

                # Link to diagnoses table
                commonPatientsDXs = df_diagnoses[(df_diagnoses.PATID.isin(patientsWithThisRX)) &
                                              (df_diagnoses.admitYear == pd.to_datetime(x.RX_START_DATE).dt.year.unique()[0])]

                return pd.Series({'diagsCrossFreq': commonPatientsDXs.DX.value_counts().index,
                                 'counts': commonPatientsDXs.DX.value_counts().values})


            meds_correlatedDiagsCoefficients = (grouped_meds.apply(patientsCrossFreqMedsDiags))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for med in meds_correlatedDiagsCoefficients.index:
                thisMedYear = meds_correlatedDiagsCoefficients.loc[med]
                thisMed = med[0]
                thisYear = med[1]
                totalCrossTab = np.sum(thisMedYear.counts)
                for diagInd in range(len(meds_correlatedDiagsCoefficients.loc[med].diagsCrossFreq.values)):
                    mytups.append((thisMedYear.diagsCrossFreq.values[diagInd], thisMedYear.counts[diagInd]/totalCrossTab))
                    multiIndex.append((thisMed, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            meds_correlatedDiagsCoefficients = (pd.DataFrame.from_records(mytups, columns=['DX','Relative_Counts'],
                                                                             index=index))

            #meds_correlatedMedsCoefficients
            def patientsCrossFreqMedsMeds(x):

                patientsWithThisRX = list(set(x.PATID.tolist()))

                # Link to labs table
                commonPatientsMeds = df_meds[(df_meds.PATID.isin(patientsWithThisRX)) &
                                              (pd.to_datetime(df_meds.RX_START_DATE).dt.year == 
                                               pd.to_datetime(x.RX_START_DATE).dt.year.unique()[0])]

                return pd.Series({'medsCrossFreq': commonPatientsMeds.JH_INGREDIENT_RXNORM_CODE.value_counts().index,
                                 'counts': commonPatientsMeds.JH_INGREDIENT_RXNORM_CODE.value_counts().values})


            meds_correlatedMedsCoefficients = (grouped_meds.apply(patientsCrossFreqMedsMeds))

            # Currently a little hacky, but seems fast

            mytups = list()
            multiIndex = list()

            for med in meds_correlatedMedsCoefficients.index:
                thisMedYear = meds_correlatedMedsCoefficients.loc[med]
                thisMed = med[0]
                thisYear = med[1]
                totalCrossTab = np.sum(thisMedYear.counts)
                for medInd in range(len(meds_correlatedMedsCoefficients.loc[med].medsCrossFreq.values)):
                    mytups.append((thisMedYear.medsCrossFreq.values[medInd], thisMedYear.counts[medInd]/totalCrossTab))
                    multiIndex.append((thisMed, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            meds_correlatedMedsCoefficients = (pd.DataFrame.from_records(mytups, columns=['JH_INGREDIENT_RXNORM_CODE','Relative_Counts'],
                                                                             index=index))

            ## MEDS TO PROCEDURES
            def patientsCrossFreqMedsProcs(x):

                patientsWithThisRX = list(set(x.PATID.tolist()))

                # Link to procs table
                commonPatientsProcs = df_procedures[df_procedures.PATID.isin(patientsWithThisRX) &
                                             (df_procedures.encounterYear == pd.to_datetime(x.RX_START_DATE).dt.year.unique()[0])]

                return pd.Series({'procsCrossFreq': commonPatientsProcs.RAW_PX.value_counts().index,
                                 'counts': commonPatientsProcs.RAW_PX.value_counts().values})

            meds_correlatedProceduresCoefficients = (grouped_meds.apply(patientsCrossFreqMedsProcs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for med in meds_correlatedProceduresCoefficients.index:
                thisMedYear = meds_correlatedProceduresCoefficients.loc[med]
                thisMed = med[0]
                thisYear = med[1]
                totalCrossTab = np.sum(thisMedYear.counts)
                for procInd in range(len(meds_correlatedProceduresCoefficients.loc[med].procsCrossFreq.values)):
                    mytups.append((thisMedYear.procsCrossFreq.values[procInd], thisMedYear.counts[procInd]/totalCrossTab))
                    multiIndex.append((thisMed, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)    
            meds_correlatedProceduresCoefficients = (pd.DataFrame.from_records(mytups, columns=['RAW_PX','Relative_Counts'],
                                                                              index=index))

            ## MEDS TO HPO
            def patientsCrossFreqMedsHPOs(x):

                patientsWithThisRX = list(set(x.PATID.tolist()))

                # Link to hpo table
                commonPatientsHPOs = df_phenotypes[(df_phenotypes.PATID.isin(patientsWithThisRX)) &
                                              (df_phenotypes.admitYear == pd.to_datetime(x.RX_START_DATE).dt.year.unique()[0])]

                return pd.Series({'hposCrossFreq': commonPatientsHPOs.HPO.value_counts().index,
                                 'counts': commonPatientsHPOs.HPO.value_counts().values})

            meds_correlatedPhenotypesCoefficients = (grouped_meds.apply(patientsCrossFreqMedsHPOs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for med in meds_correlatedPhenotypesCoefficients.index:
                thisMedYear = meds_correlatedPhenotypesCoefficients.loc[med]
                thisMed = med[0]
                thisYear = med[1]
                totalCrossTab = np.sum(thisMedYear.counts)
                for phenoInd in range(len(meds_correlatedPhenotypesCoefficients.loc[med].hposCrossFreq.values)):
                    mytups.append((thisMedYear.hposCrossFreq.values[phenoInd], thisMedYear.counts[phenoInd]/totalCrossTab))
                    multiIndex.append((thisMed, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)    
            meds_correlatedPhenotypesCoefficients = (pd.DataFrame.from_records(mytups, columns=['HPO','Relative_Counts'],
                                                                              index=index))

            return (meds_medication, meds_dosageInfo, meds_frequencyPerYear, meds_fractionOfSubjects,
                    meds_correlatedLabsCoefficients, meds_correlatedDiagsCoefficients, meds_correlatedMedsCoefficients,
                   meds_correlatedProceduresCoefficients, meds_correlatedPhenotypesCoefficients)
        
        # Make Procedures Profile
        elif profileType == 'procedures':
            procedures_code = df_procedures.RAW_PX.unique()
            procedures_count = df_procedures.RAW_PX.value_counts()

            procedures_frequencyPerYear = (df_procedures.groupby(['RAW_PX','encounterYear','PATID']).PATID.count()
                                                    .groupby(['RAW_PX','encounterYear']).mean())

            procedures_fractionOfSubjects = (np.divide(df_procedures.groupby(['RAW_PX']).PATID.nunique(),
                                            df_procedures.PATID.nunique()))

            grouped_procs = df_procedures.groupby(['RAW_PX', 'encounterYear'])

            #procs_correlatedLabsCoefficients
            def patientsAboveBelowNormalProcsLabs(x):

                patientsWithThisProc = list(set(x.PATID.tolist()))

                # Link to labs table
                abnormalPatientsLabs = df_labs[(df_labs.PATID.isin(patientsWithThisProc)) & 
                                               ((df_labs.RESULT_NUM > df_labs.range_high) | 
                                                (df_labs.RESULT_NUM < df_labs.range_low)) &
                                              (df_labs.resultYear == pd.to_datetime(x.PX_DATE).dt.year.unique()[0])]

                return pd.Series({'labsAboveBelowNorm': abnormalPatientsLabs.LAB_LOINC.value_counts().index,
                                 'counts': abnormalPatientsLabs.LAB_LOINC.value_counts().values})

            procs_correlatedLabsCoefficients = (grouped_procs.apply(patientsAboveBelowNormalProcsLabs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for proc in procs_correlatedLabsCoefficients.index:
                thisProcYear = procs_correlatedLabsCoefficients.loc[proc]
                thisProc = proc[0]
                thisYear = proc[1]
                totalCrossTab = np.sum(thisProcYear.counts)
                for labInd in range(len(procs_correlatedLabsCoefficients.loc[proc].labsAboveBelowNorm.values)):
                    mytups.append((thisProcYear.labsAboveBelowNorm.values[labInd], thisProcYear.counts[labInd]/totalCrossTab))
                    multiIndex.append((thisProc, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            procs_correlatedLabsCoefficients = (pd.DataFrame.from_records(mytups, columns=['LAB_LOINC','Relative_Counts'],
                                                                             index=index))

            #procs_correlatedDiagsCoefficients
            def patientsCrossFreqProcsDiags(x):

                patientsWithThisProc = list(set(x.PATID.tolist()))

                # Link to diagnoses table
                commonPatientsDXs = df_diagnoses[(df_diagnoses.PATID.isin(patientsWithThisProc)) &
                                              (df_diagnoses.admitYear == pd.to_datetime(x.PX_DATE).dt.year.unique()[0])]

                return pd.Series({'diagsCrossFreq': commonPatientsDXs.DX.value_counts().index,
                                 'counts': commonPatientsDXs.DX.value_counts().values})


            procs_correlatedDiagsCoefficients = (grouped_procs.apply(patientsCrossFreqProcsDiags))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for proc in procs_correlatedDiagsCoefficients.index:
                thisProcYear = procs_correlatedDiagsCoefficients.loc[proc]
                thisProc = proc[0]
                thisYear = proc[1]
                totalCrossTab = np.sum(thisProcYear.counts)
                for diagInd in range(len(procs_correlatedDiagsCoefficients.loc[proc].diagsCrossFreq.values)):
                    mytups.append((thisProcYear.diagsCrossFreq.values[diagInd], thisProcYear.counts[diagInd]/totalCrossTab))
                    multiIndex.append((thisProc, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            procs_correlatedDiagsCoefficients = (pd.DataFrame.from_records(mytups, columns=['DX','Relative_Counts'],
                                                                             index=index))

            #procs_correlatedMedsCoefficients
            def patientsCrossFreqProcsMeds(x):

                patientsWithThisProc = list(set(x.PATID.tolist()))

                # Link to labs table
                commonPatientsMeds = df_meds[(df_meds.PATID.isin(patientsWithThisProc)) &
                                              (df_meds.startYear == pd.to_datetime(x.PX_DATE).dt.year.unique()[0])]

                return pd.Series({'medsCrossFreq': commonPatientsMeds.JH_INGREDIENT_RXNORM_CODE.value_counts().index,
                                 'counts': commonPatientsMeds.JH_INGREDIENT_RXNORM_CODE.value_counts().values})


            procs_correlatedMedsCoefficients = (grouped_procs.apply(patientsCrossFreqProcsMeds))

            # Currently a little hacky, but seems fast

            mytups = list()
            multiIndex = list()

            for proc in procs_correlatedMedsCoefficients.index:
                thisProcYear = procs_correlatedMedsCoefficients.loc[proc]
                thisProc = proc[0]
                thisYear = proc[1]
                totalCrossTab = np.sum(thisProcYear.counts)
                for medInd in range(len(procs_correlatedMedsCoefficients.loc[proc].medsCrossFreq.values)):
                    mytups.append((thisProcYear.medsCrossFreq.values[medInd], thisProcYear.counts[medInd]/totalCrossTab))
                    multiIndex.append((thisProc, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            procs_correlatedMedsCoefficients = (pd.DataFrame.from_records(mytups, columns=['JH_INGREDIENT_RXNORM_CODE','Relative_Counts'],
                                                                             index=index))


            ## PROCEDURES TO PROCEDURES
            def patientsCrossFreqProcsProcs(x):

                patientsWithThisProc = list(set(x.PATID.tolist()))

                # Link to procs table
                commonPatientsProcs = df_procedures[df_procedures.PATID.isin(patientsWithThisProc) &
                                             (df_procedures.encounterYear == pd.to_datetime(x.PX_DATE).dt.year.unique()[0])]

                return pd.Series({'procsCrossFreq': commonPatientsProcs.RAW_PX.value_counts().index,
                                 'counts': commonPatientsProcs.RAW_PX.value_counts().values})

            procs_correlatedProceduresCoefficients = (grouped_procs.apply(patientsCrossFreqProcsProcs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for proc in procs_correlatedProceduresCoefficients.index:
                thisProcYear = procs_correlatedProceduresCoefficients.loc[proc]
                thisProc = proc[0]
                thisYear = proc[1]
                totalCrossTab = np.sum(thisProcYear.counts)
                for procInd in range(len(procs_correlatedProceduresCoefficients.loc[proc].procsCrossFreq.values)):
                    mytups.append((thisProcYear.procsCrossFreq.values[procInd], thisProcYear.counts[procInd]/totalCrossTab))
                    multiIndex.append((thisProc, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)    
            procs_correlatedProceduresCoefficients = (pd.DataFrame.from_records(mytups, columns=['RAW_PX','Relative_Counts'],
                                                                              index=index))

            # procedures to hpo
            def patientsCrossFreqProcsHPOs(x):

                patientsWithThisProc = list(set(x.PATID.tolist()))

                # Link to diagnoses table
                commonPatientsHPOs = df_phenotypes[(df_phenotypes.PATID.isin(patientsWithThisProc)) &
                                              (df_phenotypes.admitYear == pd.to_datetime(x.PX_DATE).dt.year.unique()[0])]

                return pd.Series({'hposCrossFreq': commonPatientsHPOs.HPO.value_counts().index,
                                 'counts': commonPatientsHPOs.HPO.value_counts().values})


            procs_correlatedPhenotypesCoefficients = (grouped_procs.apply(patientsCrossFreqProcsHPOs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for proc in procs_correlatedPhenotypesCoefficients.index:
                thisProcYear = procs_correlatedPhenotypesCoefficients.loc[proc]
                thisProc = proc[0]
                thisYear = proc[1]
                totalCrossTab = np.sum(thisProcYear.counts)
                for phenoInd in range(len(procs_correlatedPhenotypesCoefficients.loc[proc].hposCrossFreq.values)):
                    mytups.append((thisProcYear.hposCrossFreq.values[phenoInd], thisProcYear.counts[phenoInd]/totalCrossTab))
                    multiIndex.append((thisProc, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            procs_correlatedPhenotypesCoefficients = (pd.DataFrame.from_records(mytups, columns=['HPO','Relative_Counts'],
                                                                             index=index))


            return (procedures_code, procedures_count, procedures_frequencyPerYear, procedures_fractionOfSubjects,
                    procs_correlatedLabsCoefficients, procs_correlatedDiagsCoefficients, procs_correlatedMedsCoefficients,
                   procs_correlatedProceduresCoefficients, procs_correlatedPhenotypesCoefficients)
        
        # Make Diagnoses Profile
        elif profileType == 'diagnoses':
            diagnoses_code = df_diagnoses.DX.unique()
    
            diagnoses_count = df_diagnoses.DX.value_counts()

            diagnoses_frequencyPerYear = (df_diagnoses.groupby(['DX','admitYear','PATID']).PATID
                                .count().groupby(['DX','admitYear']).mean())

            diagnoses_fractionOfSubjects = (np.divide(df_diagnoses.groupby(['DX']).PATID.nunique(),
                                            df_diagnoses.PATID.nunique()))

            grouped_diags = df_diagnoses.groupby(['DX','admitYear'])

            #diags_correlatedLabsCoefficients
            def patientsAboveBelowNormalDiagsLabs(x):

                patientsWithThisDiag = list(set(x.PATID.tolist()))

                # Link to labs table
                abnormalPatientsLabs = df_labs[(df_labs.PATID.isin(patientsWithThisDiag)) & 
                                               ((df_labs.RESULT_NUM > df_labs.range_high) | 
                                                (df_labs.RESULT_NUM < df_labs.range_low)) &
                                              (df_labs.resultYear == pd.to_datetime(x.ADMIT_DATE).dt.year.unique()[0])]

                return pd.Series({'labsAboveBelowNorm': abnormalPatientsLabs.LAB_LOINC.value_counts().index,
                                 'counts': abnormalPatientsLabs.LAB_LOINC.value_counts().values})

            diags_correlatedLabsCoefficients = (grouped_diags.apply(patientsAboveBelowNormalDiagsLabs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for diag in diags_correlatedLabsCoefficients.index:
                thisDiagYear = diags_correlatedLabsCoefficients.loc[diag]
                thisDiag = diag[0]
                thisYear = diag[1]
                totalCrossTab = np.sum(thisDiagYear.counts)
                for labInd in range(len(diags_correlatedLabsCoefficients.loc[diag].labsAboveBelowNorm.values)):
                    mytups.append((thisDiagYear.labsAboveBelowNorm.values[labInd], thisDiagYear.counts[labInd]/totalCrossTab))
                    multiIndex.append((thisDiag, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            diags_correlatedLabsCoefficients = (pd.DataFrame.from_records(mytups, columns=['LAB_LOINC','Relative_Counts'],
                                                                             index=index))

            #diags_correlatedDiagsCoefficients
            def patientsCrossFreqDiagsDiags(x):

                patientsWithThisDiag = list(set(x.PATID.tolist()))

                # Link to diagnoses table
                commonPatientsDXs = df_diagnoses[(df_diagnoses.PATID.isin(patientsWithThisDiag)) &
                                              (df_diagnoses.admitYear == pd.to_datetime(x.ADMIT_DATE).dt.year.unique()[0])]

                return pd.Series({'diagsCrossFreq': commonPatientsDXs.DX.value_counts().index,
                                 'counts': commonPatientsDXs.DX.value_counts().values})


            diags_correlatedDiagsCoefficients = (grouped_diags.apply(patientsCrossFreqDiagsDiags))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for diag in diags_correlatedDiagsCoefficients.index:
                thisDiagYear = diags_correlatedDiagsCoefficients.loc[diag]
                thisDiag = diag[0]
                thisYear = diag[1]
                totalCrossTab = np.sum(thisDiagYear.counts)
                for diagInd in range(len(diags_correlatedDiagsCoefficients.loc[diag].diagsCrossFreq.values)):
                    mytups.append((thisDiagYear.diagsCrossFreq.values[diagInd], thisDiagYear.counts[diagInd]/totalCrossTab))
                    multiIndex.append((thisDiag, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            diags_correlatedDiagsCoefficients = (pd.DataFrame.from_records(mytups, columns=['DX','Relative_Counts'],
                                                                             index=index))

            #diags_correlatedMedsCoefficients
            def patientsCrossFreqDiagsMeds(x):

                patientsWithThisDiag = list(set(x.PATID.tolist()))

                # Link to labs table
                commonPatientsMeds = df_meds[(df_meds.PATID.isin(patientsWithThisDiag)) &
                                              (df_meds.startYear == pd.to_datetime(x.ADMIT_DATE).dt.year.unique()[0])]

                return pd.Series({'medsCrossFreq': commonPatientsMeds.JH_INGREDIENT_RXNORM_CODE.value_counts().index,
                                 'counts': commonPatientsMeds.JH_INGREDIENT_RXNORM_CODE.value_counts().values})


            diags_correlatedMedsCoefficients = (grouped_diags.apply(patientsCrossFreqDiagsMeds))

            # Currently a little hacky, but seems fast

            mytups = list()
            multiIndex = list()

            for diag in diags_correlatedMedsCoefficients.index:
                thisDiagYear = diags_correlatedMedsCoefficients.loc[diag]
                thisDiag = diag[0]
                thisYear = diag[1]
                totalCrossTab = np.sum(thisDiagYear.counts)
                for medInd in range(len(diags_correlatedMedsCoefficients.loc[diag].medsCrossFreq.values)):
                    mytups.append((thisDiagYear.medsCrossFreq.values[medInd], thisDiagYear.counts[medInd]/totalCrossTab))
                    multiIndex.append((thisDiag, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            diags_correlatedMedsCoefficients = (pd.DataFrame.from_records(mytups, columns=['JH_INGREDIENT_RXNORM_CODE','Relative_Counts'],
                                                                             index=index))


            ## DIAGNOSES TO PROCEDURES
            def patientsCrossFreqDiagsProcs(x):

                patientsWithThisDiag = list(set(x.PATID.tolist()))

                # Link to procs table
                commonPatientsProcs = df_procedures[df_procedures.PATID.isin(patientsWithThisDiag) &
                                             (df_procedures.encounterYear == pd.to_datetime(x.ADMIT_DATE).dt.year.unique()[0])]

                return pd.Series({'procsCrossFreq': commonPatientsProcs.RAW_PX.value_counts().index,
                                 'counts': commonPatientsProcs.RAW_PX.value_counts().values})

            diags_correlatedProceduresCoefficients = (grouped_diags.apply(patientsCrossFreqDiagsProcs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for diag in diags_correlatedProceduresCoefficients.index:
                thisDiagYear = diags_correlatedProceduresCoefficients.loc[diag]
                thisDiag = diag[0]
                thisYear = diag[1]
                totalCrossTab = np.sum(thisDiagYear.counts)
                for procInd in range(len(diags_correlatedProceduresCoefficients.loc[diag].procsCrossFreq.values)):
                    mytups.append((thisDiagYear.procsCrossFreq.values[procInd], thisDiagYear.counts[procInd]/totalCrossTab))
                    multiIndex.append((thisDiag, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)    
            diags_correlatedProceduresCoefficients = (pd.DataFrame.from_records(mytups, columns=['RAW_PX','Relative_Counts'],
                                                                              index=index))

            #diags_correlatedPhenotypesCoefficients
            def patientsCrossFreqDiagsHPOs(x):

                patientsWithThisDiag = list(set(x.PATID.tolist()))

                # Link to diagnoses table
                commonPatientsHPOs = df_phenotypes[(df_phenotypes.PATID.isin(patientsWithThisDiag)) &
                                              (df_phenotypes.admitYear == pd.to_datetime(x.ADMIT_DATE).dt.year.unique()[0])]

                return pd.Series({'hposCrossFreq': commonPatientsHPOs.HPO.value_counts().index,
                                 'counts': commonPatientsHPOs.HPO.value_counts().values})


            diags_correlatedPhenotypesCoefficients = (grouped_diags.apply(patientsCrossFreqDiagsHPOs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for diag in diags_correlatedPhenotypesCoefficients.index:
                thisDiagYear = diags_correlatedPhenotypesCoefficients.loc[diag]
                thisDiag = diag[0]
                thisYear = diag[1]
                totalCrossTab = np.sum(thisDiagYear.counts)
                for phenoInd in range(len(diags_correlatedPhenotypesCoefficients.loc[diag].hposCrossFreq.values)):
                    mytups.append((thisDiagYear.hposCrossFreq.values[phenoInd], thisDiagYear.counts[phenoInd]/totalCrossTab))
                    multiIndex.append((thisDiag, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            diags_correlatedPhenotypesCoefficients = (pd.DataFrame.from_records(mytups, columns=['HPO','Relative_Counts'],
                                                                             index=index))

            return (diagnoses_code, diagnoses_count, diagnoses_frequencyPerYear, diagnoses_fractionOfSubjects,
                    diags_correlatedLabsCoefficients, diags_correlatedDiagsCoefficients, diags_correlatedMedsCoefficients,
                   diags_correlatedProceduresCoefficients, diags_correlatedPhenotypesCoefficients)
        
        # Make Phenotypes Profile
        elif profileType == 'phenotypes':
            phenotypes_code = df_phenotypes.HPO.unique()

            phenotypes_count = df_phenotypes.HPO.value_counts()

            phenotypes_frequencyPerYear = (df_phenotypes.groupby(['HPO','admitYear','PATID']).PATID
                                .count().groupby(['HPO','admitYear']).mean())

            phenotypes_fractionOfSubjects = (np.divide(df_phenotypes.groupby(['HPO']).PATID.nunique(),
                                            df_phenotypes.PATID.nunique()))

            grouped_phenotypes = df_phenotypes.groupby(['HPO','admitYear'])

            #diags_correlatedLabsCoefficients
            def patientsAboveBelowNormalDiagsLabs(x):

                patientsWithThisDiag = list(set(x.PATID.tolist()))

                # Link to labs table
                abnormalPatientsLabs = df_labs[(df_labs.PATID.isin(patientsWithThisDiag)) & 
                                               ((df_labs.RESULT_NUM > df_labs.range_high) | 
                                                (df_labs.RESULT_NUM < df_labs.range_low)) &
                                              (df_labs.resultYear == pd.to_datetime(x.ADMIT_DATE).dt.year.unique()[0])]

                return pd.Series({'labsAboveBelowNorm': abnormalPatientsLabs.LAB_LOINC.value_counts().index,
                                 'counts': abnormalPatientsLabs.LAB_LOINC.value_counts().values})

            phenos_correlatedLabsCoefficients = (grouped_phenotypes.apply(patientsAboveBelowNormalDiagsLabs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for diag in phenos_correlatedLabsCoefficients.index:
                thisDiagYear = phenos_correlatedLabsCoefficients.loc[diag]
                thisDiag = diag[0]
                thisYear = diag[1]
                totalCrossTab = np.sum(thisDiagYear.counts)
                for labInd in range(len(phenos_correlatedLabsCoefficients.loc[diag].labsAboveBelowNorm.values)):
                    mytups.append((thisDiagYear.labsAboveBelowNorm.values[labInd], thisDiagYear.counts[labInd]/totalCrossTab))
                    multiIndex.append((thisDiag, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            phenos_correlatedLabsCoefficients = (pd.DataFrame.from_records(mytups, columns=['LAB_LOINC','Relative_Counts'],
                                                                             index=index))

            #diags_correlatedDiagsCoefficients
            def patientsCrossFreqDiagsDiags(x):

                patientsWithThisDiag = list(set(x.PATID.tolist()))

                # Link to diagnoses table
                commonPatientsDXs = df_diagnoses[(df_diagnoses.PATID.isin(patientsWithThisDiag)) &
                                              (df_diagnoses.admitYear == pd.to_datetime(x.ADMIT_DATE).dt.year.unique()[0])]

                return pd.Series({'diagsCrossFreq': commonPatientsDXs.DX.value_counts().index,
                                 'counts': commonPatientsDXs.DX.value_counts().values})


            phenos_correlatedDiagsCoefficients = (grouped_phenotypes.apply(patientsCrossFreqDiagsDiags))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for diag in phenos_correlatedDiagsCoefficients.index:
                thisDiagYear = phenos_correlatedDiagsCoefficients.loc[diag]
                thisDiag = diag[0]
                thisYear = diag[1]
                totalCrossTab = np.sum(thisDiagYear.counts)
                for diagInd in range(len(phenos_correlatedDiagsCoefficients.loc[diag].diagsCrossFreq.values)):
                    mytups.append((thisDiagYear.diagsCrossFreq.values[diagInd], thisDiagYear.counts[diagInd]/totalCrossTab))
                    multiIndex.append((thisDiag, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            phenos_correlatedDiagsCoefficients = (pd.DataFrame.from_records(mytups, columns=['DX','Relative_Counts'],
                                                                             index=index))

            #diags_correlatedMedsCoefficients
            def patientsCrossFreqDiagsMeds(x):

                patientsWithThisDiag = list(set(x.PATID.tolist()))

                # Link to labs table
                commonPatientsMeds = df_meds[(df_meds.PATID.isin(patientsWithThisDiag)) &
                                              (df_meds.startYear == pd.to_datetime(x.ADMIT_DATE).dt.year.unique()[0])]

                return pd.Series({'medsCrossFreq': commonPatientsMeds.JH_INGREDIENT_RXNORM_CODE.value_counts().index,
                                 'counts': commonPatientsMeds.JH_INGREDIENT_RXNORM_CODE.value_counts().values})


            phenos_correlatedMedsCoefficients = (grouped_phenotypes.apply(patientsCrossFreqDiagsMeds))

            # Currently a little hacky, but seems fast

            mytups = list()
            multiIndex = list()

            for diag in phenos_correlatedMedsCoefficients.index:
                thisDiagYear = phenos_correlatedMedsCoefficients.loc[diag]
                thisDiag = diag[0]
                thisYear = diag[1]
                totalCrossTab = np.sum(thisDiagYear.counts)
                for medInd in range(len(phenos_correlatedMedsCoefficients.loc[diag].medsCrossFreq.values)):
                    mytups.append((thisDiagYear.medsCrossFreq.values[medInd], thisDiagYear.counts[medInd]/totalCrossTab))
                    multiIndex.append((thisDiag, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            phenos_correlatedMedsCoefficients = (pd.DataFrame.from_records(mytups, columns=['JH_INGREDIENT_RXNORM_CODE','Relative_Counts'],
                                                                             index=index))


            ## DIAGNOSES TO PROCEDURES
            def patientsCrossFreqDiagsProcs(x):

                patientsWithThisDiag = list(set(x.PATID.tolist()))

                # Link to procs table
                commonPatientsProcs = df_procedures[df_procedures.PATID.isin(patientsWithThisDiag) &
                                             (df_procedures.encounterYear == pd.to_datetime(x.ADMIT_DATE).dt.year.unique()[0])]

                return pd.Series({'procsCrossFreq': commonPatientsProcs.RAW_PX.value_counts().index,
                                 'counts': commonPatientsProcs.RAW_PX.value_counts().values})

            phenos_correlatedProceduresCoefficients = (grouped_phenotypes.apply(patientsCrossFreqDiagsProcs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for diag in phenos_correlatedProceduresCoefficients.index:
                thisDiagYear = phenos_correlatedProceduresCoefficients.loc[diag]
                thisDiag = diag[0]
                thisYear = diag[1]
                totalCrossTab = np.sum(thisDiagYear.counts)
                for procInd in range(len(phenos_correlatedProceduresCoefficients.loc[diag].procsCrossFreq.values)):
                    mytups.append((thisDiagYear.procsCrossFreq.values[procInd], thisDiagYear.counts[procInd]/totalCrossTab))
                    multiIndex.append((thisDiag, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)    
            phenos_correlatedProceduresCoefficients = (pd.DataFrame.from_records(mytups, columns=['RAW_PX','Relative_Counts'],
                                                                              index=index))

            #diags_correlatedPhenotypesCoefficients
            def patientsCrossFreqDiagsHPOs(x):

                patientsWithThisDiag = list(set(x.PATID.tolist()))

                # Link to diagnoses table
                commonPatientsHPOs = df_phenotypes[(df_phenotypes.PATID.isin(patientsWithThisDiag)) &
                                              (df_phenotypes.admitYear == pd.to_datetime(x.ADMIT_DATE).dt.year.unique()[0])]

                return pd.Series({'hposCrossFreq': commonPatientsHPOs.HPO.value_counts().index,
                                 'counts': commonPatientsHPOs.HPO.value_counts().values})


            phenos_correlatedPhenotypesCoefficients = (grouped_phenotypes.apply(patientsCrossFreqDiagsHPOs))

            # Currently a little hacky, but seems fast
            mytups = list()
            multiIndex = list()

            for diag in phenos_correlatedPhenotypesCoefficients.index:
                thisDiagYear = phenos_correlatedPhenotypesCoefficients.loc[diag]
                thisDiag = diag[0]
                thisYear = diag[1]
                totalCrossTab = np.sum(thisDiagYear.counts)
                for phenoInd in range(len(phenos_correlatedPhenotypesCoefficients.loc[diag].hposCrossFreq.values)):
                    mytups.append((thisDiagYear.hposCrossFreq.values[phenoInd], thisDiagYear.counts[phenoInd]/totalCrossTab))
                    multiIndex.append((thisDiag, thisYear))

            index = pd.MultiIndex.from_tuples(multiIndex)

            phenos_correlatedPhenotypesCoefficients = (pd.DataFrame.from_records(mytups, columns=['HPO','Relative_Counts'],
                                                                             index=index))

            return (phenotypes_code, phenotypes_count, phenotypes_frequencyPerYear, phenotypes_fractionOfSubjects,
                    phenos_correlatedLabsCoefficients, phenos_correlatedDiagsCoefficients, phenos_correlatedMedsCoefficients,
                   phenos_correlatedProceduresCoefficients, phenos_correlatedPhenotypesCoefficients)
        else:
            print("Please provide profileType as 'labs', 'medications', 'procedures', 'diagnoses', or 'phenotypes'")
            
    except ValueError:
        print("Please provide profileType as 'labs', 'medications', 'procedures', 'diagnoses', or 'phenotypes'")