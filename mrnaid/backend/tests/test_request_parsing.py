from unittest import main, TestCase
from unittest.mock import MagicMock as Mock
from common.utils.RequestParser import RequestParser
from json import dumps


class MockRequest:
    def __init__(self,
        config,
        dinucleotides, uridine_depletion,
        faster_MFE_algorithm, file_name, sequences):
        self.config = config
        self.dinucleotides = dinucleotides
        self.uridine_depletion = uridine_depletion
        self.faster_MFE_algorithm = faster_MFE_algorithm
        self.file_name = file_name
        self.sequences = sequences
        

    def __repr__(self):
        return f"MockRequest({self.__dict__})"

    def __str__(self):
        return f"MockRequest({self.__dict__})"

    def __call__(self, *args, **kwargs):
        return MockRequest(*args, **kwargs)
    
    def get_data(self):
        return dumps({
            "config": self.config,
            "dinucleotides": self.dinucleotides,
            "uridine_depletion": self.uridine_depletion,
            "faster_MFE_algorithm": self.faster_MFE_algorithm,
            "file_name": self.file_name,
            "sequences": self.sequences,
            "optimization_criterion": "codon_usage",
        })

# Configure MockRequest to have the attributes and methods you need
MockRequest.return_value.get_data.return_value = dumps({
    "config": {
        "avoided_motifs": [
            "EcoRI"
        ],
        "codon_usage_frequency_threshold": 0.1,
        "max_GC_content": 0.7,
        "min_GC_content": 0.3,
        "GC_window_size": 100,
        "organism": "h_sapiens",
        "number_of_sequences": 3,
        "entropy_window": 30,
        "optimizationAlgorithm": "nil",
        "importedCodonTable": "undefined"
    },
    "uridine_depletion": True,
    "optimization_criterion": "codon_usage",
    "precise_MFE_algorithm": True,
    "file_name": "results",
    "sequences": {
        "five_end_flanking_sequence": "AGAACUAAG",
        "gene_of_interest": "CGCUUUGCCAAAGUUGGCCAGAAGCCCACGCAUUAAUUCUAGUGUGUGGGGAUGUGAAGAAUGUCCCGGAAAAGAUGUCCAUUGCGCUGGAGACCCUGGAUUGGCUACGCAUGAUCUCCUUCCCUAGACUUCGGUUCCUUCUAAGAGGUAUACUCCACAGCCCAGUUGGCCCAUAGAGACGCCUACUAACGAAACGGUCACUAUAGACCAAAAACUACGAUGUCCGCGAUGGAGAUACGCUGACUUACUACCCAUUCAUUCCGCGACAGUAGAGAAAAGACGACGGACUAUCAUUCUGAAGAGUCGUUCAAGUAGGAUUCCGGCCGACCUGUACCUAAAUGGCGUCCGUCCGGGUGUCCUUGUAACUGCCUAGAAACAGGAGCAGGACGGGAAGUCCUCGUAGUCCUAUCGUGGUUCUGGCCUCAAAUCCGCCCCGAAGAGUGAGCAACGACUCGUCGUACAUUUAGUCCACGGGAGAUCUGAACAAUCUCGGUGGUCUACCUUGCUUGUUUCGACACUUUCUUCCAUCCGAGCAAAUAUCGGGCUACCCACUUGAGGUAUGAGGGACCC",
        "three_end_flanking_sequence": "UUGAGUGCAAUUUCCAAUGCAGCCUAC"
    }
})


class TestRequestParser(TestCase):
    
    def test_correct_parsing(self):
        """Test if the request JSON is parsed correctly"""
        config = {
            "avoid_codons": ["ATG"],
            "avoided_motifs": [],
            "codon_usage_frequency_threshold": 0.1,
            "max_GC_content": 0.7,
            "GC_window_size": 30,
            "min_GC_content": 0.3,
            "organism": 'h_sapiens',
            "entropy_window": 6,
            "number_of_sequences": 5
        }
        dinucleotides = False
        uridine_depletion = False
        faster_MFE_algorithm = False
        filename = 'test'
        sequences = {"five_end_flanking_sequence": "ACG",
                    "gene_of_interest": "AGCAAAGAAGGCGGA",
                    "three_end_flanking_sequence": "ACG"}
        request = MockRequest(config=config,
                            dinucleotides=dinucleotides, uridine_depletion=uridine_depletion,
                            faster_MFE_algorithm=faster_MFE_algorithm, file_name=filename, sequences=sequences)
        
        parser = RequestParser(request)
        params = parser.parse()
        print(params.input_mRNA)
        # try:
        #     params = parser.parse()
        # except Exception as e:
        #     print(f'Exception is {e}')
        #     success = False
        # else:
        #     success = True

        # assert success
        assert params.input_mRNA == "AGCAAAGAAGGCGGA"
        assert len(params.avoid_motifs) == 0
        assert len(params.avoid_codons) == 1 and params.avoid_codons[0] == "ATG"
        assert params.max_GC_content == 0.7
        assert params.min_GC_content == 0.3
        assert params.GC_window_size == 30
        assert params.usage_threshold == 0.1
        assert not params.uridine_depletion
        assert params.organism == 'h_sapiens'
        assert params.entropy_window == 6
        assert params.number_of_sequences == 5
        assert params.filename == 'test'
        assert params.mfe_method == 'RNAfold'
        assert not params.dinucleotides

if __name__ == '__main__':
    main()