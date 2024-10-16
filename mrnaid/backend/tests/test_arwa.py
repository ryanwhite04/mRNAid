from unittest import main, TestCase
from common.arw_mrna.src.awalk import adaptive_random_walk, WalkConfig
from common.arw_mrna.src.protein import CodonFrequencyTable, random_cds
from common.arw_mrna.src.objective_functions import make_cai_threshold_obj, make_cai_and_efe_obj, make_cai_and_aup_obj, CAIThresholdObjectiveConfig

class TestARWA(TestCase):
        
    def test_adaptive_random_walk(self):
        aa_seq = "MVSKGEELFTGVVPIL"
        freq_table = CodonFrequencyTable("homosapiens.txt")
        objective = make_cai_and_efe_obj(
            CAIThresholdObjectiveConfig(
                freq_table=freq_table,
                cai_threshold=0.8,
                cai_exp_scale=1.0,
                verbose=False
            )
        )
        cds = ['AUG', 'GUU', 'UCC', 'AAG', 'GGA', 'GAA', 'GAA', 'CUC', 'UUC', 'ACU', 'GGU', 'GUU', 'GUU', 'CCG', 'AUC', 'UUG', 'GUU', 'GAA', 'CUC', 'GAU', 'GGA', 'GAU', 'GUA', 'AAU', 'GGA', 'CAC']
        
        walk_config = WalkConfig(
            aa_seq,
            freq_table,
            objective,
            steps=1,
            init_cds=cds,
            seed=42
        )
        generator = adaptive_random_walk(walk_config)
        initial = next(generator)
        self.assertEqual(initial.step, 0)
        self.assertAlmostEqual(initial.best_fitness, 16.45672606, places=6)
    
if __name__ == '__main__':
    main()