from unittest import main, TestCase

# from common.arw_mrna.src.protein import CodonFrequencyTable
from common.arw_mrna.src.protein import CodonFrequencyTable

class TestCodonTable(TestCase):
    
    def test_get_homosapien(self):
        table = CodonFrequencyTable("common/arw_mrna/codon_tables/homosapiens.txt")
        self.assertEqual(table.get_codon_freq("ACA"), 0.28) # 614523
        self.assertEqual(table.get_aa_max_freq("T"), 0.36) # 768147
        self.assertEqual(table.get_codons("T"), {"ACA", "ACC", "ACG", "ACU"})
        self.assertEqual(table.get_aa("ACA"), "T")
        self.assertEqual(table.max_codons(), 6)
        self.assertEqual(round(table.codon_adaption_weight("ACA"), 2), 0.78) # 0.8
        self.assertEqual(round(table.codon_adaption_weight("ACC"), 2), 1)   # 1
        self.assertEqual(round(table.codon_adaption_weight("ACG"), 2), 0.31) # 0.32
        self.assertEqual(round(table.codon_adaption_weight("ACU"), 2), 0.69) # 0.69
        self.assertEqual(round(table.codon_adaptation_index(["ACA", "ACC", "ACG", "ACU"]), 2), 0.64) # 0.65
        self.assertEqual(round(table.log_codon_adaptation_index(["ACA", "ACC", "ACG", "ACU"]), 2), -0.45) # -0.43

if __name__ == '__main__':
    main()