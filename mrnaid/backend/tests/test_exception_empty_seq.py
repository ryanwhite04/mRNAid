from unittest import main, TestCase
from tests.test_request_parsing import MockRequest
from common.utils.Exceptions import EmptySequenceError
from common.utils.RequestParser import RequestParser


class TestRequestParser(TestCase):

    def setUp(self):
        self.example_parameters = {
            "config": {
                "avoid_codons": [],
                "avoided_motifs": [],
                "codon_usage_frequency_threshold": 0.1,
                "max_GC_content": 0.6,
                "GC_window_size": 30,
                "min_GC_content": 0.4,
                "organism": 'h_sapiens',
                "entropy_window": 6,
                "number_of_sequences": 5
            },
            "dinucleotides": False,
            "uridine_depletion": True,
            "faster_MFE_algorithm": False,
            "file_name": "test"
        }

    def test_empty_seq(self):
        sequences = {
            "five_end_flanking_sequence": "ACG",
            "gene_of_interest": "",
            "three_end_flanking_sequence": "ACG"
        }
        request = MockRequest(**self.example_parameters, sequences=sequences)
        parser = RequestParser(request)
        with self.assertRaises(EmptySequenceError):
            parser.parse()

    def test_empty_five_end(self):
        sequences = {
            "five_end_flanking_sequence": "",
            "gene_of_interest": "AAAGGGCCCTTT",
            "three_end_flanking_sequence": "ACG"
        }
        request = MockRequest(**self.example_parameters, sequences=sequences)
        parser = RequestParser(request)
        with self.assertRaises(EmptySequenceError):
            parser.parse()

    def test_empty_three_end(self):
        sequences = {
            "five_end_flanking_sequence": "ATG",
            "gene_of_interest": "AAAGGGCCCTTT",
            "three_end_flanking_sequence": ""
        }
        request = MockRequest(**self.example_parameters, sequences=sequences)
        parser = RequestParser(request)
        with self.assertRaises(EmptySequenceError):
            parser.parse()


if __name__ == '__main__':
    main()
