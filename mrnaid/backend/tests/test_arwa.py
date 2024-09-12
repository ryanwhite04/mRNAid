from unittest import main, TestCase
from common.arw_mrna.src.awalk import adaptive_random_walk, WalkConfig
from common.arw_mrna.src.protein import CodonFrequencyTable
from common.arw_mrna.src.objective_functions import make_cai_threshold_obj, make_cai_and_efe_obj, make_cai_and_aup_obj, CAIThresholdObjectiveConfig

class TestARWA(TestCase):
    
    def test_adaptive_random_walk(self):
        walk_config = WalkConfig(
            aa_seq="MVSKGEELFTGVVPILVELDGDVNGH",
            freq_table=CodonFrequencyTable("homosapiens.txt"),
            objective=make_cai_and_efe_obj,
            steps=100,
            verbose=True,
            init_cds=None
        )
        result = adaptive_random_walk(walk_config)
    
if __name__ == '__main__':
    main()