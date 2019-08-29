def getSubdemographicsTables(user, passwd, cohort='All', schema='dbo', medianEncounterYear=2019, sex='All', race='All', age_low='All', 
                             age_high=None ):
    """Extract data from PCORnet database and clean to prepare for Clinical Profile calculation
    
    Keyword arguments:
    user, psswd -- for access to PCORnet data 
    cohort -- the name of the table / view from PCORnet on Azure (default 'All', meaning no disease cohort)
    schema -- schema that the table lives in on Azure (default 'dbo')
    medianEncounterYear -- used to calculate the age of people in the database, based on DOB (default 2019)
    sex -- sex to extract as part of the demographic profile (default 'All', options: 'All','M','F')
    race -- race to extract as part of the demographic profile (default 'All', options: 'All', 'White or Caucasian', 
    'Black or African American', 'Other')
    age_low -- low end of the age range to extract as part of the demographic profile (default 'All', 
    meaning specify no age restriction)
    age_high -- high end of the age range to extract as part of the demographic profile (default None, meaning specify no
    age restriction)
    
    Returns dataframes for use in calculateAnyProfile function
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
    
    driver='FreeTDS'
    tds_ver='8.0'

    host_ip='esmpmdbdev6.database.windows.net'
    db_port=1433
    db='ClinicalProfile'

    conn_str=('DRIVER={};SERVER={};PORT={};DATABASE={};UID={};PWD={};TDS_VERSION={}'.format(
    driver, host_ip, db_port, db, user, passwd, tds_ver))

    engine = sqlalchemy.create_engine('mssql+pyodbc:///?odbc_connect='+urllib.parse.quote(conn_str))
    
    if (cohort == 'All'):
        df_sub_demographics = pd.read_sql_table('DEMOGRAPHIC',str(engine.url), index_col='PATID', schema=schema)
    else:
        query = "select * from [{0}].[{1}]".format(schema, cohort)
        df_sub_demographics = pd.read_sql_query(query, engine)
    
    df_sub_demographics['race_code'] = df_sub_demographics.RACE.map({'01':'Other','02':'Other',
                                                                     '03':'Black or African American',
                                                            '04':'Other','05':'White or Caucasian','06':'Other',
                                                            '07':'Other','NI':'Other','UN':'Other','OT':'Other'})
    if (sex != 'All'):
        # grab gender
        df_sub_demographics = df_sub_demographics[df_sub_demographics.SEX == sex]
        
    if (race != 'All'):
        # grab race
        df_sub_demographics = df_sub_demographics[df_sub_demographics.race_code == race]
        
    if (age_low != 'All'):
        # grab age
        dob_ub = medianEncounterYear - float(age_low)
        dob_lb = medianEncounterYear - float(age_high)
        df_sub_demographics = (df_sub_demographics[
            (pd.to_datetime(df_sub_demographics.BIRTH_DATE).dt.year >= dob_lb) & 
            (pd.to_datetime(df_sub_demographics.BIRTH_DATE).dt.year <= dob_ub)]) 
        
    # Labs
    query = """select demo.PATID, ENCOUNTERID, LAB_LOINC, RESULT_DATE, RESULT_NUM, LOINC_SHORTNAME, LOINC_UNIT,
                    RANGE_LOW, RANGE_HIGH from [{0}].[{1}] as lab
                    inner join [{2}].[{3}] as demo
                    on lab.PATID = demo.PATID""".format('dbo', 'vw_pc_labs', schema, cohort)
    df_labs = pd.read_sql_query(query, engine)

    query = ("""select * from [dbo].[jh_loinc] """)
    df_loincinfo = pd.read_sql_query(query, engine)

    df_labs = df_sub_demographics.merge(df_labs, how='left', left_on='PATID', right_on='PATID')
    df_labs_full = df_labs.merge(df_loincinfo, how='left', left_on='LAB_LOINC', right_on='Loinc_Code')

    df_labs_full['resultYear'] = pd.to_datetime(df_labs_full.RESULT_DATE).dt.year

    df_labs_full['range_high'] = (df_labs_full.RANGE_HIGH.mask(df_labs_full.RANGE_HIGH.eq('None')).dropna()
                                 .astype('str').str.replace(',','').replace('', np.nan).replace(' ', np.nan).astype('float'))
    df_labs_full['range_low'] = (df_labs_full.RANGE_LOW.mask(df_labs_full.RANGE_LOW.eq('None')).dropna()
                                 .astype('str').str.replace(',','').replace('', np.nan).replace(' ', np.nan).astype('float'))
    df_labs_full = df_labs_full.drop(['RANGE_HIGH','RANGE_LOW'],axis=1)
#     labInfo = df_labs_full[['RESULT_NUM','LAB_LOINC','range_high','range_low','resultYear','PATID', 'LOINC_SHORTNAME']]

    
    # Meds
    query = """select demo.PATID, RX_START_DATE, JH_INGREDIENT_RXNORM_CODE, RX_DOSE_ORDERED, RX_QUANTITY, RX_ROUTE
                    from [{0}].[{1}] as med
                    inner join [{2}].[{3}] as demo
                    on med.PATID = demo.PATID""".format('dbo', 'PRESCRIBING', schema, cohort)
    df_meds = pd.read_sql_query(query, engine)
    
    df_meds_full = df_sub_demographics.merge(df_meds, how='left', left_on='PATID', right_on='PATID')
    
    df_meds_full['startYear'] = pd.to_datetime(df_meds_full.RX_START_DATE).dt.year
#     rxInfo = df_meds_full[['JH_INGREDIENT_RXNORM_CODE', 'PATID', 'startYear']]
    
    # Procedures
    query = """select demo.PATID, ENCOUNTERID, RAW_PX, PX_DATE
                    from [{0}].[{1}] as p
                    inner join [{2}].[{3}] as demo
                    on p.PATID = demo.PATID""".format('dbo', 'PROCEDURES', schema, cohort)
    df_procedures = pd.read_sql_query(query, engine)
    
    df_procedures_full = df_sub_demographics.merge(df_procedures, how='left', left_on='PATID', right_on='PATID')
    
    df_procedures_full['encounterYear'] = pd.to_datetime(df_procedures_full.PX_DATE).dt.year
#     procInfo = df_procedures_full[['RAW_PX','PATID', 'encounterYear']]
    
    # Diagnoses
    query = """select demo.PATID, ENCOUNTERID, DX, ADMIT_DATE
                    from [{0}].[{1}] as diag
                    inner join [{2}].[{3}] as demo
                    on diag.PATID = demo.PATID
                    where diag.DX_TYPE = '10'""".format('dbo', 'DIAGNOSIS', schema, cohort)
    df_diagnoses = pd.read_sql_query(query, engine)
    df_diagnoses_full = df_sub_demographics.merge(df_diagnoses, how='left', left_on='PATID', right_on='PATID')
    
    df_diagnoses_full['admitYear'] = pd.to_datetime(df_diagnoses_full.ADMIT_DATE).dt.year
#     diagInfo = df_diagnoses_full[['DX','PATID','admitYear']]
    
    # HPO
    hpoMapping = pd.read_csv('HPO_ICD10_nostring.txt', sep=' ', header=None)
    hpoMapping.drop(2,axis=1,inplace=True)
    hpoMapping.columns = ['ICD10','HPO']
    
    df_phenotypes_full = df_diagnoses_full.merge(hpoMapping, left_on='DX', right_on='ICD10', how='inner')
    df_phenotypes_full['admitYear'] = pd.to_datetime(df_phenotypes_full.ADMIT_DATE).dt.year
#     phenoInfo = df_phenotypes_full[['HPO','PATID', 'admitYear']]
        
    return (df_labs_full, df_meds_full, df_procedures_full, df_diagnoses_full, df_phenotypes_full)