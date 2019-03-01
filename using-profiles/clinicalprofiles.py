import requests
import pandas as pd


DEFAULT_KEYS = [
    'resourceType', 'id', 'meta', 'text',
    'identifier', 'status', 'population',
    'cohort', 'date', 'source', 'reporter'
]


class ClinicalProfileServer:

    def __init__(self, base_url):
        self.base_url = base_url

        self._data = None
        self._populate_data()

    def _populate_data(self):
        res = requests.get(f"{self.base_url}/baseR4/ClinicalProfile?_pretty=true")
        if res.ok:
            self._data = res.json()
        else:
            raise ValueError(
                f"Failed to fetch from {self.base_url}/baseR4/ClinicalProfile?_pretty=true"
            )

    def keys(self):
        return [
            self._data['entry'][i]['resource']['id']
            for i in range(len(self._data['entry']))
        ]

    def profiles(self):
        return self.keys()

    def __getitem__(self, key):
        for entry in self._data['entry']:
            if entry['resource']['id'] == key:
                return ClinicalProfile.from_dict(entry)



class Variable:

    def __init__(self, name, data):
        self._name = name
        self._data = data

    @property
    def name(self):
        return self._name['display']

    @property
    def code(self):
        return self._name['code']

    @property
    def system(self):
        return self._name['system']

    @property
    def fraction_of_subjects(self):
        return self._data.get('fractionOfSubjects', None)

    @property
    def count(self):
        return self._data.get('count', None)

    def __repr__(self):
        return f"<Variable '{self.name}'>"

    @staticmethod
    def from_dict(data: dict) -> "Variable":
        return Variable(data['code'][0]['coding'][0], data)

class PhenotypeVariable(Variable):

    @property
    def phenotypes(self):
        return {
            x['code']['coding'][0]['display']: x.get('coefficient', None)
            for x in self._data['correlatedPhenotypes']['entry']
        }

    @property
    def phenotype_codes(self):
        return {
            x['code']['coding'][0]['name']: x['code']['coding'][0]['display']
            for x in self._data['correlatedPhenotypes']['entry']
        }

    @property
    def raw_phenotype(self):
        return self._data['correlatedPhenotypes']

    def __len__(self):
        return self.phenotypes

    def __getitem__(self, key):
        return self.phenotypes[key]

    def keys(self):
        return self.phenotypes.keys()

    @staticmethod
    def from_dict(data: dict) -> "PhenotypeVariable":
        return PhenotypeVariable(data['code'][0]['coding'][0], data)


class LabVariable(Variable):
    # Statistics:
    @property
    def distribution(self):
        return self._data['scalarDistribution']
    @property
    def units(self):
        return self._data['scalarDistribution']['units']
    @property
    def min(self):
        return self._data['scalarDistribution']['min']
    @property
    def max(self):
        return self._data['scalarDistribution']['max']
    @property
    def mean(self):
        return self._data['scalarDistribution']['mean']
    @property
    def stdDev(self):
        return self._data['scalarDistribution']['stdDev']
    @property
    def fraction_above_normal(self):
        return self._data['scalarDistribution']['fractionAboveNormal']
    @property
    def fraction_below_normal(self):
        return self._data['scalarDistribution']['fractionBelowNormal']

    @staticmethod
    def from_dict(data: dict) -> "LabVariable":
        return LabVariable(data['code'][0]['coding'][0], data)


class MedicationVariable(Variable):

    @staticmethod
    def from_dict(data: dict) -> "MedicationVariable":
        return MedicationVariable(data['dosage']['route'][0]['coding'][0], data)

class DiagnosisVariable(Variable):
    pass

class ProcedureVariable(Variable):
    pass

class ClinicalProfile:

    def __init__(self, name, data):
        self.name = name
        self._data = data

        self._resource = self._data['resource']

        self.last_updated = pd.Timestamp(self._resource['meta']['lastUpdated']).to_pydatetime()
        self.date = pd.Timestamp(self._resource['date']).to_pydatetime()
        self.version = self._resource['meta']['versionId']
        self.url = self._data['fullUrl']
        self.status = self._resource['status']

        self.population = self._resource['population']
        self.cohort = self._resource['cohort']
        self.source = self._resource['source']
        self.reporter = self._resource['reporter']

        self.phenotypes = self._resource["hpo"]
        self.labs = self._resource["lab"]
        self.medications = self._resource["medication"]
        self.diagnoses = self._resource["diagnosis"]
        self.procedures = self._resource["procedure"]


    def get_phenotype_variables(self) -> "List[Variable]":
        return {
            PhenotypeVariable.from_dict(x).name: PhenotypeVariable.from_dict(x)
            for x in self.phenotypes
        }

    def get_phenotype_codes(self) -> "List[Variable]":
        return {
            PhenotypeVariable.from_dict(x).code: PhenotypeVariable.from_dict(x)
            for x in self.phenotypes
        }


    def get_lab_variables(self) -> "List[Variable]":
        return {
            LabVariable.from_dict(x).name: LabVariable.from_dict(x)
            for x in self.labs
        }

    def get_lab_codes(self) -> "List[Variable]":
        return {
            LabVariable.from_dict(x).code: LabVariable.from_dict(x)
            for x in self.labs
        }


    def get_medication_variables(self) -> "List[Variable]":
        return {
            MedicationVariable.from_dict(x).name: MedicationVariable.from_dict(x)
            for x in self.medications
        }

    def get_medication_codes(self) -> "List[Variable]":
        return {
            MedicationVariable.from_dict(x).code: MedicationVariable.from_dict(x)
            for x in self.medications
        }


    def get_diagnosis_variables(self) -> "List[Variable]":
        return {
            DiagnosisVariable.from_dict(x).name: DiagnosisVariable.from_dict(x)
            for x in self.diagnoses
        }

    def get_diagnosis_codes(self) -> "List[Variable]":
        return {
            DiagnosisVariable.from_dict(x).code: DiagnosisVariable.from_dict(x)
            for x in self.diagnoses
        }


    def get_procedure_variables(self) -> "List[Variable]":
        return {
            ProcedureVariable.from_dict(x).name: ProcedureVariable.from_dict(x)
            for x in self.procedures
        }

    def get_procedure_codes(self) -> "List[Variable]":
        return {
            ProcedureVariable.from_dict(x).code: ProcedureVariable.from_dict(x)
            for x in self.procedures
        }


    def __repr__(self):
        return f"<ClinicalProfile '{self.name}'>"

    @staticmethod
    def from_dict(data: dict) -> 'ClinicalProfile':
        return ClinicalProfile(data['resource']['id'], data)
