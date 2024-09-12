# common/Evaluation.py
```python
import math
from typing import List, Dict
from typing import Tuple

import RNA
from Bio.Seq import Seq
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex, SynonymousCodons
from billiard import Pool
from dnachisel import Location
from dnachisel import MatchTargetCodonUsage
from python_codon_tables import get_codons_table

from common.objectives.Codon_pair_usage import MatchTargetCodonPairUsage
from common.objectives.Dinucleotide_usage import MatchTargetPairUsage
from common.objectives.MFE import MFE
from common.utils.Datatypes import OptimizationParameters
from common.utils.Datatypes import SequenceProperties, EvaluationResults
from common.utils.Logger import MyLogger

# setting up a logger
logger = MyLogger(__name__)


class Evaluation(object):
    """ This class allows evaluate the final results by calculating sequence parameters and ranking the final results"""

    def __init__(self, optimized_seqs: List[str], parameters: OptimizationParameters):
        self.input_seq = parameters.input_mRNA.upper()
        self.optimized_seqs = optimized_seqs
        self.max_GC_content = parameters.max_GC_content
        self.min_GC_content = parameters.min_GC_content
        self.GC_window_size = parameters.GC_window_size
        self.usage_threshold = parameters.usage_threshold
        self.uridine_depletion = parameters.uridine_depletion
        self.organism = parameters.organism
        self.location = parameters.location
        self.entropy_window = parameters.entropy_window
        self.five_end = parameters.five_end.replace('T', 'U')
        self.three_end = parameters.three_end.replace('T', 'U')
        self.evaluation_dict = {}
        self.dinucleotides = parameters.optimization_criterion == 'dinucleotides'
        self.codon_pair = parameters.optimization_criterion == 'codon_pair'

    def score(self, sequence: str, gc_ratio: float, mfe5: float, mfe_total: float) -> float:
        """
        Use sequence properties to calculate a score for ranking.
        :param sequence: RNA sequence
        :param gc_ratio: GC ratio in the sequence
        :param mfe5: MFE of the 5' end (calculated by RNAfold)
        :param mfe_total: MFE of the full sequence
        :return: ranking score
        """

        # MFE at 5' end
        mfe_5_score = math.exp(mfe5 / 100)

        # total MFE score
        mfe_score = math.exp(mfe_total / 1000)

        # Score reflecting how far we are from the target codon usage
        codon_score = 1 - (self.get_codon_comparison_organism(sequence) / (len(sequence) / 3))
        # minimum: 0
        # maximum: 1 (when normalized by length in codons) (maximum absolute value)

        # count each U at third position in codon, normalize to codon number
        u_depletion = sum(1 for i in range(0, len(sequence) - 2, 3) if sequence[i + 2] in ["U", "T"]) / (
                len(sequence) / 3)
        # minimum: 0 (no codons have U at third pos)
        # maximum: 1 (all codons have U)

        # GC content score
        gc_score = (gc_ratio - self.min_GC_content) / (self.max_GC_content - self.min_GC_content)

        # Dinucleotide frequency score
        # If a custom table is used in Objectives.py, it needs to be used here, too.
        dinucleotides = MatchTargetPairUsage()
        pair_frequencies = dinucleotides.calculate_freqs(sequence)
        di_score = self.differences_from_table(dinucleotides.pair_usage_table, pair_frequencies)
        dinucleotide_score = 1 - di_score  # This way the score should be between 0 and 1, where 1 is the best
        # normalized by number of all dinucleotides

        # Codon pair frequency score
        codon_pairs = MatchTargetCodonPairUsage()
        codon_pair_freqs = codon_pairs.calculate_freqs(sequence)
        pair_score = self.differences_from_table(codon_pairs.pair_usage_table, codon_pair_freqs)
        codon_pair_score = 1 - pair_score

        # CAI score (coincides with CAI value)
        cai_score = self.CAI(sequence)

        # Score for ranking the output
        # Uridine depletion, GC content, CAI, MFE global, MFE 5' end scores are considered
        w1 = 1  # weight for mfe5 score
        w2 = 8  # weight for gc_sore
        w3 = 10 if self.uridine_depletion else 0  # weight for Uridine depletion score
        w4 = 0 if not (
                    self.codon_pair and self.dinucleotides) else 0  # not used now, can be included by weight modification
        w5 = 0 if self.dinucleotides else 0  # not used now, can be included by weight modification
        w6 = 0 if self.codon_pair else 0  # not used now, can be included by weight modification
        w7 = 7  # weight for CAI score
        w8 = 3  # weight for total MFE score
        score = (w1 * mfe_5_score + w2 * gc_score + w3 * u_depletion + w4 * codon_score +
                 w5 * dinucleotide_score + w6 * codon_pair_score + w7 * cai_score +
                 w8 * mfe_score) / (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8)
        return score

    def differences_from_table(self, table: Dict[str, float], frequencies: Dict[str, float]) -> float:
        score = 0
        for pair, freq in frequencies.items():
            score_for_pair = abs(freq - frequencies[pair])
            score += score_for_pair
        return score / len(table)

    def get_codon_dict(self, sequence: str) -> Dict[str, dict]:
        """Create a comparison between codon frequencies in the sequence and in the chosen organism."""
        sequence = sequence.replace('U', 'T')
        comparison = MatchTargetCodonUsage(species=self.organism, location=Location(0, len(sequence), 1))
        codons = self.get_codons(sequence)

        codons_positions, aa_comparisons = comparison.compare_frequencies(codons)
        for codon in codons:
            aa = Seq(codon).translate()
            aa_comparisons[aa][codon]["count"] = aa_comparisons[aa][codon].get("count", 0) + 1
        """
        Returns dictionary in format:
         {
           "K": {
             "total": 6,
             "AAA": {
                 "sequence": 1.0,
                 "table": 0.7,
                 "count": 6
             },
             ...
           },
           "D": ...
         }
        """
        return aa_comparisons

    def get_codon_comparison_organism(self, sequence: str) -> float:
        """
        For each codon, calculate absolute difference between frequency in the sequence and in the organism, weighted by
        total number of given codon in the sequence. Return sum of these weighted differences.
        """
        aa_comparisons = self.get_codon_dict(sequence)
        codon_score = 0
        for aa, data in aa_comparisons.items():
            total = data.pop("total")
            for codon, codon_freq in data.items():
                frequency_diff = codon_freq["sequence"] - codon_freq["table"]
                count = aa_comparisons[aa][codon].get("count", 0)
                codon_score += count * abs(frequency_diff)
        return codon_score

    def get_codons(self, sequence: str) -> List[str]:
        """Return a list of codons in the sequence."""
        return [
            sequence[3 * i: 3 * (i + 1)]
            for i in range(int(len(sequence) / 3))
        ]

    def _get_max_frequency(self, codon, usage_table):
        """
        Get frequency of the most frequent synonymous codon to codon for that amino-acid - the scaling factor in the
        RCSU value (CAI is the geometric mean of RCSU values of codons)
        """
        synonymous_codons = [codons for (aa, codons) in SynonymousCodons.items() if codon in codons][0]
        max_freq = max([usage_table[codon] for codon in synonymous_codons])
        return max_freq

    def codon_table_to_biopython_index(self, codon_table: dict) -> dict:
        """Transforms table from python-codon-table library to biopython index format"""
        amino_acids = ['*', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                       'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

        dicts = [codon_table[aa] for aa in amino_acids]
        simple_usage_table = {key: val for dic in dicts for (key, val) in dic.items()}
        return {key: val / self._get_max_frequency(key, simple_usage_table) for key, val in simple_usage_table.items()}

    def CAI(self, sequence: str) -> float:
        """Calculates Codon Adaptation Index for a given sequence"""
        sequence = sequence.replace('U', 'T')
        cai_estimator = CodonAdaptationIndex()
        codon_table = get_codons_table(self.organism)
        index = self.codon_table_to_biopython_index(codon_table)
        cai_estimator.set_cai_index(index)
        return cai_estimator.cai_for_gene(sequence)

    def get_seq_properties(self, tag_seq: Tuple) -> SequenceProperties:
        """
        Calculate properties for the input sequence, return them as a dictionary.
        :param tag_seq: tuple of tag (tag of the sequence), seq (RNA sequence)
        :return: Dictionary of parameters
        """
        # The entire sequence consisting of fixed and variable parts
        full_seq = self.five_end + tag_seq[1] + self.three_end
        rna_sequence = full_seq.replace('T', 'U')
        dna_sequence = full_seq.replace('U', 'T')
        base_ratio = dict()
        seq_len = len(full_seq)
        # Count each nucleotide type
        a_s = sum(full_seq[x] in ['A', 'a'] for x in range(seq_len))
        ts = sum(full_seq[x] in ['U', 'u', 'T', 't'] for x in range(seq_len))
        gs = sum(full_seq[x] in ['G', 'g'] for x in range(seq_len))
        cs = sum(full_seq[x] in ['C', 'c'] for x in range(seq_len))
        # Calculate their ratio
        base_ratio['A'] = a_s / seq_len
        base_ratio['TU'] = ts / seq_len
        base_ratio['G'] = gs / seq_len
        base_ratio['C'] = cs / seq_len
        gc_ratio = base_ratio['G'] + base_ratio['C']
        at_ratio = 1 - gc_ratio

        mfe_estimator = MFE()

        rna_structure, mfe = RNA.fold(full_seq)
        # 5' MFE is calculated up to this index
        index_for_mfe_evaluation = len(self.five_end) + self.entropy_window
        _, mfe5 = RNA.fold(full_seq[:index_for_mfe_evaluation])
        _, mfe_total = RNA.fold(full_seq[index_for_mfe_evaluation:])
        score = self.score(tag_seq[1], gc_ratio, mfe5, mfe_total)

        # calculating CAI
        cai = self.CAI(tag_seq[1])

        logger.info("----------" + str(tag_seq[0]))
        for base, ratio in base_ratio.items():
            logger.info(base + "\t" + str(ratio))
        logger.info(f">>AT Ratio: {at_ratio}")
        logger.info(f">>GC Ratio: {gc_ratio}")
        logger.info(f">>MFE at 5': {mfe5}")
        logger.info(f">>Total MFE: {mfe}")
        logger.info(f">>CAI: {cai}")
        logger.info(f"Score: {score}")

        properties = SequenceProperties(seqID=tag_seq[0],
                                        RNASeq=rna_sequence,
                                        DNASeq=dna_sequence,
                                        RNA_structure=rna_structure,
                                        A_ratio=base_ratio['A'],
                                        TU_ratio=base_ratio['TU'],
                                        G_ratio=base_ratio['G'],
                                        C_ratio=base_ratio['C'],
                                        AT_ratio=at_ratio,
                                        GC_ratio=gc_ratio,
                                        MFE=mfe,
                                        MFE_5=mfe5,
                                        score=score,
                                        CAI=cai)

        return properties

    def get_evaluation(self) -> Dict:
        """
        For the set of sequences, create a dictionary containing basic properties for each of them.
        :return: Dictionary of sequence properties
        """
        tag_seqs = [('opt_' + str(i), seq) for i, seq in enumerate(self.optimized_seqs)]

        pool = Pool()
        optimized_results = pool.map(self.get_seq_properties, tag_seqs)

        optimized_results_sorted = sorted(optimized_results, key=lambda k: k.score, reverse=True)

        evaluation_results = EvaluationResults(input=self.get_seq_properties(('input', self.input_seq)),
                                               optimized=optimized_results_sorted)

        return evaluation_results.dict()

```
    
# common/OptimizationProblems.py
```python
from common.constraints.Constraints import Constraints
from dnachisel import DnaOptimizationProblem
from dnachisel.Location import Location
from common.objectives.Objectives import Objectives
from common.utils.Datatypes import OptimizationParameters
from common.utils.Exceptions import SpeciesError


def initialize_optimization_problem(parameters: OptimizationParameters) -> DnaOptimizationProblem:
    """
    Combines all constraints and objectives together to create a list of DNA Optimization problems.
    :param parameters: optimization parameters
    :return: optimization problem object
    """

    constraints_creator = Constraints(parameters.avoid_motifs, parameters.uridine_depletion,
                                      parameters.min_GC_content, parameters.max_GC_content, parameters.organism,
                                      parameters.usage_threshold,
                                      Location(0, len(parameters.input_DNA)))
    try:
        constraints = constraints_creator.create_constraints()
    except FileNotFoundError as e:
        raise SpeciesError("No data provided for the organism {}".format(parameters.organism))

    # Create objectives
    objectives_creator = Objectives(parameters.entropy_window, parameters.organism,
                                    Location(0, len(parameters.input_DNA)), parameters.min_GC_content,
                                    parameters.max_GC_content, parameters.GC_window_size, parameters.five_end,
                                    parameters.mfe_method, parameters.optimization_criterion)

    try:
        objectives = objectives_creator.create_objectives()
    except FileNotFoundError as e:
        raise SpeciesError("No data provided for the organism {}".format(parameters.organism))

    return DnaOptimizationProblem(sequence=parameters.input_DNA, constraints=constraints, objectives=objectives)

```
    
# common/OptimizationTask.py
```python
from dnachisel import DnaOptimizationProblem
from common.utils.Logger import MyLogger
from common.utils.Exceptions import OptimizationFailedError

logger = MyLogger(__name__)


def optimization_task(problem: DnaOptimizationProblem) -> str:
    """
    Try to resolve constraints and optimize given problem, catch exceptions, if any.
    :param problem: DnaOptimizationProblem with specified constraints and objectives
    :return: optimized sequence
    """

    try:
        # Resolve constrains
        problem.resolve_constraints()
        # Optimize
        problem.optimize()

        print(problem.constraints_text_summary())
        print(problem.objectives_text_summary())

        # Optimized sequence
        optimized = problem.sequence.replace('T', 'U')
        print("Optimized mRNA sequence:{}".format(optimized))
        logger.debug("Optimized mRNA sequence:{}".format(optimized))
    except Exception as e:
        logger.error(e.message)
        raise OptimizationFailedError()
    else:
        return optimized

```
    
# common/arw_mrna/src/awalk.py
```python
""" Implements adaptive random walks on CDSs. """
from typing import Callable, List, Optional
from .protein import CodonFrequencyTable, random_cds
import random
import dataclasses


@dataclasses.dataclass
class WalkConfig:
    aa_seq: str
    freq_table: CodonFrequencyTable
    objective: Callable[[List[str]], float]
    steps: int
    init_cds: Optional[List[str]] = None
    verbose: bool = False


@dataclasses.dataclass
class WalkResult:
    cds: List[str]
    fitness: float


def adaptive_random_walk(config: WalkConfig):
    """Runs ARW to maximize objective function on CDSs."""
    cds = random_cds(
        config.aa_seq, config.freq_table) if config.init_cds is None else config.init_cds
    prev_fitness = config.objective(cds)
    if config.verbose:
        print(f"Initial CDS: {cds}")
    for step in range(config.steps):
        mutidx = random.randint(0, len(cds) - 1)
        while len(config.freq_table.get_codons(config.aa_seq[mutidx])) == 1:
            mutidx = random.randint(0, len(cds) - 1)
        mutcodons = list(config.freq_table.get_codons(
            config.aa_seq[mutidx])-set([cds[mutidx]]))
        mutcodon = random.choice(mutcodons)
        new_cds = cds[:mutidx] + [mutcodon] + cds[mutidx+1:]
        new_fitness = config.objective(new_cds)
        if new_fitness > prev_fitness:
            cds = new_cds
            prev_fitness = new_fitness
            if config.verbose:
                print(f"New CDS: {new_cds}")
        if config.verbose:
            print(
                f"Step: {step}, Fitness: {new_fitness}, Best Fitness: {prev_fitness}")
    return WalkResult(cds, prev_fitness)


def adaptive_random_walk_generator(config: WalkConfig):
    """Runs ARW to maximize objective function on CDSs."""
    cds = random_cds(
        config.aa_seq, config.freq_table) if config.init_cds is None else config.init_cds
    prev_fitness, measures = config.objective(cds)
    if config.verbose:
        yield {
            "type":"initial",
            "cds": cds,
            "fitness": prev_fitness,
            "measures": measures
        }
    for step in range(config.steps):
        mutidx = random.randint(0, len(cds) - 1)
        while len(config.freq_table.get_codons(config.aa_seq[mutidx])) == 1:
            mutidx = random.randint(0, len(cds) - 1)
        mutcodons = list(config.freq_table.get_codons(
            config.aa_seq[mutidx]) - set([cds[mutidx]]))
        mutcodon = random.choice(mutcodons)
        new_cds = cds[:mutidx] + [mutcodon] + cds[mutidx+1:]
        new_fitness, measures = config.objective(new_cds)
        if new_fitness > prev_fitness:
            cds = new_cds
            prev_fitness = new_fitness
            if config.verbose:
                yield {
                    "type": "new_cds",
                    "cds": new_cds,
                    "fitness": new_fitness,
                    "step": step,
                    "best_fitness": prev_fitness,
                    "measures": measures
                }
        if config.verbose:
            yield {
                "type": "progress",
                "step": step,
                "fitness": new_fitness,
                "best_fitness": prev_fitness,
                "measures": measures
            }
    yield {
        "type": "final",
        "cds": cds,
        "step": config.steps - 1,
        "fitness": prev_fitness,
        "measures": measures
    }
```
    
# common/arw_mrna/src/objective_functions.py
```python
"""Implementations of objective functions for use in ARWs."""
from typing import Callable, List
import math
import dataclasses
from .vienna import cds_to_rna, average_unpaired, ensemble_free_energy
from .protein import CodonFrequencyTable


@dataclasses.dataclass
class CAIThresholdObjectiveConfig:
    """Configuration for CAI threshold objective function."""
    freq_table: CodonFrequencyTable
    cai_threshold: float = 0.8
    cai_exp_scale: float = 1.0
    verbose: bool = False


def make_cai_threshold_obj(config: CAIThresholdObjectiveConfig) -> Callable[[List[str]], float]:
    """Optimises CAI up to the threshold"""
    def obj(cds: List[str]) -> float:
        cai = config.freq_table.codon_adaptation_index(cds)
        cai_penalty = math.exp(
            max(0, config.cai_threshold-cai)*config.cai_exp_scale)-1
        if config.verbose:
            print(f"-- Obj fn log. CAI: {cai}")
        fitness = -cai_penalty
        measures = {"CAI": cai}
        return fitness, measures
    return obj


def make_cai_and_aup_obj(config: CAIThresholdObjectiveConfig) -> Callable[[List[str]], float]:
    """Optimises CAI and AUP: (1-aup)-(e^(max(0,threshold-cai)*scale)-1)"""
    def obj(cds: List[str]) -> float:
        rna_seq = cds_to_rna(cds)
        cai = config.freq_table.codon_adaptation_index(cds)
        cai_penalty = math.exp(
            max(0, config.cai_threshold-cai)*config.cai_exp_scale)-1
        aup = average_unpaired(rna_seq)
        ensemble_paired_prob = 1-aup
        if config.verbose:
            print(f"-- Obj fn log. CAI: {cai}, AUP: {aup}")
        fitness = ensemble_paired_prob - cai_penalty
        measures = {"CAI": cai, "AUP": aup}
        return fitness, measures
    return obj


def make_cai_and_efe_obj(config: CAIThresholdObjectiveConfig) -> Callable[[List[str]], float]:
    """Optimises CAI and EFE: -efe*(1/e^(max(0,threshold-cai)*scale))"""
    def obj(cds: List[str]) -> float:
        rna_seq = cds_to_rna(cds)
        cai = config.freq_table.codon_adaptation_index(cds)
        cai_penalty = math.exp(
            max(0, config.cai_threshold-cai)*config.cai_exp_scale)
        efe = ensemble_free_energy(rna_seq)
        if config.verbose:
            print(f"-- Obj fn log. CAI: {cai}, EFE: {efe}")
        fitness = -efe*(1/cai_penalty)
        measures = {"CAI": cai, "EFE": efe}
        return fitness, measures
    return obj

```
    
# common/arw_mrna/src/protein.py
```python
""" This module contains code for dealing with amino acid sequence and coding sequence (CDS) data. """
import math
import dataclasses

AA_SINGLE_LETTER = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "End": "*",  # Stop codon
}


def is_valid_aa_letter(aa):
    return aa in AA_SINGLE_LETTER.values() and aa != '*'


class CodonFrequencyTable:
    def __init__(self, table_path):
        file = open(table_path, 'r')
        self.codons = set()
        self.codon_to_aa = {}
        self.aa_to_codons = {}
        self.codon_freq = {}
        self.aa_max_freq = {}
        # The format uses T instead of U
        self.nt_map = {'A': 'A', 'T': 'U', 'C': 'C', 'G': 'G'}
        for line in file:
            tokens = line.strip(" \n").split()
            if len(tokens) < 3:
                continue
            aa = tokens[0]
            aa = AA_SINGLE_LETTER[aa]
            codon = ''.join([self.nt_map[nt] for nt in tokens[1]])
            freq = round(float(tokens[2]))
            self.codons.add(codon)
            self.codon_to_aa[codon] = aa
            if aa not in self.aa_to_codons:
                self.aa_to_codons[aa] = set()
            self.aa_to_codons[aa].add(codon)
            self.codon_freq[codon] = freq
            if aa not in self.aa_max_freq:
                self.aa_max_freq[aa] = 0
            self.aa_max_freq[aa] = max(self.aa_max_freq[aa], freq)
        file.close()

    def get_codon_freq(self, codon):
        return self.codon_freq[codon]

    def get_aa_max_freq(self, aa):
        return self.aa_max_freq[aa]

    def get_codons(self, aa) -> set[str]:
        return self.aa_to_codons[aa]

    # Maximum number of codons for a single amino acid
    def max_codons(self) -> int:
        return max(len(self.get_codons(aa)) for aa in self.aa_to_codons)

    def get_aa(self, codon):
        return self.codon_to_aa[codon]

    def codon_adaption_weight(self, codon):
        return self.get_codon_freq(codon) / self.get_aa_max_freq(self.get_aa(codon))

    def codon_adaptation_index(self, cds) -> float:
        cai = 1
        for codon in cds:
            cai *= self.codon_adaption_weight(codon)
        return cai**(1/len(cds))

    def log_codon_adaptation_index(self, cds):
        cai = 0
        for codon in cds:
            cai += math.log(self.codon_adaption_weight(codon))
        return cai / len(cds)


@dataclasses.dataclass
class UniProtAASeq:
    seq: str
    uniprot_name: str
    protein_name: str


def read_cds(tsv_path):
    cds = []
    first = True
    for line in open(tsv_path, 'r'):
        if first:
            first = False
            continue
        tokens = line.strip(" \n").split("\t")
        if len(tokens) != 7:
            assert False, "Invalid CDS file format"
        if any(not is_valid_aa_letter(aa) for aa in tokens[6]):
            continue
        cds.append(UniProtAASeq(tokens[6], tokens[1], tokens[2]))
    return cds


def random_cds(aa_seq, freq_table) -> list[str]:
    import random
    cds = []
    for aa in aa_seq:
        cds.append(random.choice(list(freq_table.aa_to_codons[aa])))
    return cds


def encode_cds_one_hot(cds: list[str], freq_table: CodonFrequencyTable) -> list[list[float]]:
    width = freq_table.max_codons()
    res = []
    for codon in cds:
        enc = [0.0 for _ in range(width)]
        codons = sorted(freq_table.get_codons(freq_table.get_aa(codon)))
        enc[codons.index(codon)] = 1.0
        res.append(enc)
    return res


def main():
    """Test the codon frequency table class."""
    path = "../data/codon_tables/homosapiens.txt"
    ct = CodonFrequencyTable(path)
    print(ct.get_codon_freq("AUG"))

    cds = ["AUG", "GGC", "UUG", "UCC", "CGG", "AGC", "GAG"]
    cai = ct.codon_adaptation_index(cds)
    print(''.join(cds), cai)

    cds = ["AUG", "CCC", "CCC", "GGG", "GGC", "UUA", "AAA"]
    cai = ct.codon_adaptation_index(cds)
    print(''.join(cds), cai)

    cds = ["AUG", "GUU", "AAA", "GUA", "GUA", "GGG", "GUA"]
    cai = ct.codon_adaptation_index(cds)
    print(''.join(cds), cai)


if __name__ == "__main__":
    main()

```
    
# common/arw_mrna/src/run_awalk.py
```python
"""A script to run an adaptive walk to find an mRNA coding sequence (CDS) for a given protein."""
import pickle
import protein
import awalk
import objective_functions as objectives
import argparse
import vienna


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--steps', type=int, default=1000,
                    help='Number of steps in the adaptive walk')
    ap.add_argument('--aa_seq', type=str, required=True,
                    help='Amino acid sequence to find a CDS for. A string using the standard one-letter code.')
    ap.add_argument('--stability', type=str, choices=['aup', 'efe', 'none'], default='aup',
                    help='Stability objective to use. Set to none to only use CAI')
    ap.add_argument('--freq_table_path', type=str,
                    default='../codon_tables/homosapiens.txt', help='Path to codon frequency table')
    ap.add_argument('--cai_threshold', type=float, default=0.8,
                    help='Objective function forces CAI to be at least this')
    ap.add_argument('--cai_exp_scale', type=float, default=1.0,
                    help='Scaling factor for CAI. Increase to make CAI more important')
    ap.add_argument("--save_path", type=str, default=None,
                    help="The path to save the result. Saved in pickle format.")
    ap.add_argument("--verbose", action="store_true", help="Log all progress")
    ap.add_argument("--load_path", type=str, default=None,
                    help="Loads the initial CDS from the given pickle file. If not specified, the initial CDS is randomly generated.")
    args = ap.parse_args()

    # Load frequency table
    freq_table = protein.CodonFrequencyTable(args.freq_table_path)

    obj_config = objectives.CAIThresholdObjectiveConfig(
        freq_table, args.cai_threshold, args.cai_exp_scale, verbose=args.verbose)

    # Get obj function
    if args.stability == 'aup':
        obj = objectives.make_cai_and_aup_obj(obj_config)
    elif args.stability == 'efe':
        obj = objectives.make_cai_and_efe_obj(obj_config)
    elif args.stability == 'none':
        obj = objectives.make_cai_threshold_obj(obj_config)

    # Load initial CDS if specified
    init_cds = None
    if args.load_path is not None:
        # Read previous result as pickle
        with open(args.load_path, "rb") as f:
            init_cds = pickle.load(f).cds

    # Create walk config
    walk_config = awalk.WalkConfig(
        args.aa_seq, freq_table, obj, args.steps, init_cds=init_cds, verbose=args.verbose)

    # Run walker
    res = awalk.adapative_random_walk(walk_config)

    # Output results
    cai = freq_table.codon_adaptation_index(res.cds)
    aup, efe = vienna.aup_and_efe(vienna.cds_to_rna(res.cds))
    print(
        f"Adaptive walk result for {args.aa_seq} with {args.stability} objective and {args.steps} steps:")
    print(
        f"Result CDS: {res.cds} \n\t Fitness: {res.fitness}, CAI: {cai}, AUP: {aup}, EFE: {efe}")

    # Save the result as pickle
    if args.save_path is not None:
        with open(args.save_path, "wb") as f:
            pickle.dump(res, f)


if __name__ == '__main__':
    main()

```
    
# common/arw_mrna/src/vienna.py
```python
"""Interface to ViennaRNA package"""
from typing import Tuple
import RNA


class ViennaContext:
    def __init__(self, rna, temp=None, dangles=2, noLPs=False) -> None:
        md = RNA.md()
        md.uniq_ML = 1
        md.dangles = dangles
        md.noLP = noLPs
        if temp is not None:
            md.temperature = temp
        self.fc = RNA.fold_compound(rna, md)
        self.pf_computed = False

    def __ensure_pf(self):
        if self.pf_computed:
            return
        _, mfe_energy = self.fc.mfe()
        self.fc.exp_params_rescale(mfe_energy)
        _, efe = self.fc.pf()
        self.efe = efe
        self.pf_computed = True

    def free_energy(self, ss):
        return self.fc.eval_structure(ss)

    def prob(self, ss):
        self.__ensure_pf()
        return self.fc.pr_structure(ss)

    def make_bppt(self) -> list[list[float]]:
        self.__ensure_pf()
        bpp = self.fc.bpp()
        sz = self.fc.length
        res = [[0.0 for _ in range(sz)] for _ in range(sz)]
        for i in range(sz):
            for j in range(sz):
                if j < i:
                    res[i][j] = bpp[j+1][i+1]
                elif i < j:
                    res[i][j] = bpp[i+1][j+1]
        for i in range(sz):
            res[i][i] = 1-sum(res[i])
        return res

    def ensemble_free_energy(self):
        self.__ensure_pf()
        return self.efe

    def subopt(self, energy_delta):
        sub = self.fc.subopt(int(energy_delta*100), sorted=0)
        return [s.structure for s in sub]

    def ensemble_defect(self, ss):
        self.__ensure_pf()
        return self.fc.ensemble_defect(ss)

    def mfe(self):
        return self.fc.mfe()[0]

    def psample(self, samples=1, redundant=True):
        self.__ensure_pf()
        return self.fc.pbacktrack(samples, RNA.PBACKTRACK_DEFAULT if redundant else RNA.PBACKTRACK_NON_REDUNDANT)


def ensemble_unpaired(bppt: list[list[float]]):
    return sum([bppt[i][i] for i in range(len(bppt))])


def average_unpaired(rna_seq: str) -> float:
    ctx = ViennaContext(rna_seq)
    return ensemble_unpaired(ctx.make_bppt())/len(rna_seq)


def ensemble_free_energy(rna_seq: str) -> float:
    ctx = ViennaContext(rna_seq)
    return ctx.ensemble_free_energy()

def aup_and_efe(rna_seq: str) -> Tuple[float, float]:
    ctx = ViennaContext(rna_seq)
    return ensemble_unpaired(ctx.make_bppt())/len(rna_seq), ctx.ensemble_free_energy()


def cds_to_rna(cds):
    return ''.join(cds)

```
    
# common/constraints/Constraints.py
```python
from typing import List, Union, Tuple

from dnachisel import Specification, AvoidPattern, EnzymeSitePattern, EnforceGCContent, EnforceTranslation, \
    AvoidRareCodons
from common.constraints.UridineDepletion import UridineDepletion


class Constraints(object):

    def __init__(self, avoid_motifs: List[str], uridine_depletion: bool,
                 min_GC_content: int, max_GC_content: int, organism: str, usage_threshold: float,
                 location: Union[None, Tuple]) -> None:
        self.avoid_motifs = avoid_motifs
        self.uridine_depletion = uridine_depletion
        self.min_GC_content = min_GC_content
        self.max_GC_content = max_GC_content
        self.organism = organism
        self.usage_threshold = usage_threshold
        self.location = location

    def create_constraints(self) -> List[Specification]:
        """
        Create a list of constraints based on parameters.
        EnforceTranslation is added always, the others only if the corresponding parameter is set.
        """
        constraints_list = []
        for motif in self.avoid_motifs:
            try:
                constraints_list.append(AvoidPattern(EnzymeSitePattern(motif), location=self.location))
            except KeyError:
                motif = motif.replace('U', 'T')
                constraints_list.append(AvoidPattern(motif, location=self.location))
        constraints_list.append(EnforceGCContent(mini=self.min_GC_content, maxi=self.max_GC_content, location=self.location))
        if self.uridine_depletion:
            constraints_list.append(UridineDepletion(location=self.location))
        constraints_list.append(EnforceTranslation(genetic_table='Standard',
                                                   location=self.location))  # need different table for different species
        if self.usage_threshold:
            constraints_list.append(
                AvoidRareCodons(min_frequency=self.usage_threshold, species=self.organism, location=self.location))

        return constraints_list

```
    
# common/constraints/UridineDepletion.py
```python
from dnachisel import DnaOptimizationProblem
from dnachisel.Location import Location
from dnachisel.Specification import SpecEvaluation
from dnachisel.builtin_specifications.CodonSpecification import CodonSpecification


class UridineDepletion(CodonSpecification):
    """
    Avoid codon with U in third position in the sequence.
    """

    def __init__(self, location=None, boost=1.0):
        self.boost = boost
        self.location = Location.from_data(location)

    def initialized_on_problem(self, problem: DnaOptimizationProblem, role="constraint"):
        """Get translation from the sequence if it is not already set."""
        return self._copy_with_full_span_if_no_location(problem)

    def localized_on_window(self, new_location, start_codon, end_codon):
        return self.copy_with_changes(location=new_location)

    def evaluate(self, problem: DnaOptimizationProblem) -> SpecEvaluation:
        location = (
            self.location
            if self.location is not None
            else Location(0, len(problem.sequence))
        )
        subsequence = location.extract_sequence(problem.sequence)
        errors_locations = [
            Location(i + location.start, i + 3 + location.start, 1)  # Safe to assume strand is always 1?
            for i in range(0, len(subsequence) - 2, 3)
            if subsequence[i + 2] == "T"
        ]
        return SpecEvaluation(
            self,
            problem,
            score=-len(errors_locations),
            locations=errors_locations,
            message="All OK."
            if len(errors_locations) == 0
            else "U in third position found at indices %s" % errors_locations,
        )

    def __str__(self) -> str:
        """Represent."""
        return "UridineDepletion"

```
    
# common/constraints/__init__.py
```python

```
    
# common/objectives/Codon_pair_usage.py
```python
import json
import os
from statistics import mean
from typing import Dict, List

from dnachisel import (DnaOptimizationProblem, SpecEvaluation, Location, Specification)
from common.utils.Logger import MyLogger

logger = MyLogger(__name__)

DEFAULT_BACKEND_OBJECTIVES_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
BACKEND_OBJECTIVES_DATA = os.environ.get('BACKEND_OBJECTIVES_DATA', DEFAULT_BACKEND_OBJECTIVES_DATA)


class MatchTargetCodonPairUsage(Specification):
    """ This class is used to specify the objective to match the target
    usage of codon pairs in organism. The pair usage table is taken
    from CoCoPUTs database"""

    def __init__(self, pair_usage_table: Dict[str, float] = None, location: Location = None, species: str = 'h_sapiens'):
        super().__init__()

        if pair_usage_table is None:
            with open(os.path.join(BACKEND_OBJECTIVES_DATA, f'codon_pair_usage_{species}.json'),
                      'r') as myfile:
                data = myfile.read()
            self.pair_usage_table = json.loads(data)
        else:
            self.pair_usage_table = pair_usage_table

        self.location = Location.from_data(location)

    def evaluate(self, problem: DnaOptimizationProblem) -> SpecEvaluation:
        """Return discrepancy score for the problem"""
        sequence = problem.sequence
        calculated_freqs = self.calculate_freqs(sequence)
        score = self.calculate_score(calculated_freqs)
        breaches = self.get_overrepresented_codon_pairs(sequence, calculated_freqs)

        return SpecEvaluation(
            self, problem,
            score=-score,
            locations=breaches,
            message=f"Codon pair score: {score}"
        )

    def __str__(self) -> str:
        """String representation."""
        return "MatchTargetCodonPairUsage"

    def calculate_freqs(self, sequence: str) -> Dict[str, float]:
        """
        Calculate a dictionary of frequencies for each possible codon pair
        """
        number_of_codon_pairs = len(sequence)/3 - 1
        calculated_freqs = {}

        counts = {}
        for i in range(0, len(sequence) - 3, 3):
            counts[sequence[i:i + 6]] = counts.get(sequence[i:i + 6], 0) + 1
        for pair, _ in counts.items():
            calculated_freqs[pair] = counts.get(pair, 0) / number_of_codon_pairs

        return calculated_freqs

    def calculate_score(self, calculated_freqs: Dict[str, float]) -> float:
        ratios = []
        for pair, fr in calculated_freqs.items():
            table_frequency = self.pair_usage_table[pair]
            ratios.append(table_frequency / fr)
        return mean([abs(1 - ratio) for ratio in ratios])

    def get_overrepresented_codon_pairs(self, sequence: str, calculated_freqs: Dict[str, float]) -> List[Location]:
        pairs = [sequence[i:i + 6] for i in range(0, len(sequence) - 3, 3)]
        pair_positions = {pair: [] for pair in pairs}
        for i in range(0, len(sequence) - 3, 3):
            pair = sequence[i:i + 6]
            pair_positions[pair].append(Location(i, i+6, strand=1))
        nonoptimal = []
        for pair, freq in calculated_freqs.items():
            if freq > self.pair_usage_table[pair]:
                nonoptimal += pair_positions[pair]
        return nonoptimal

    def localized(self, location, problem=None):
        return self.copy_with_changes(location=location)


```
    
# common/objectives/Dinucleotide_usage.py
```python
import json
import os
from statistics import mean
from typing import Dict, List

from dnachisel import (DnaOptimizationProblem, SpecEvaluation, Location, Specification)
from dnachisel.biotools import group_nearby_indices
from common.utils.Logger import MyLogger

logger = MyLogger(__name__)

DEFAULT_BACKEND_OBJECTIVES_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
BACKEND_OBJECTIVES_DATA = os.environ.get('BACKEND_OBJECTIVES_DATA', DEFAULT_BACKEND_OBJECTIVES_DATA)

class MatchTargetPairUsage(Specification):
    """ This class is used to specify the objective to match the target
    usage of nucleotide pairs in organism. The pair usage table is taken
    from CoCoPUTs database: https://bigd.big.ac.cn/databasecommons/database/id/6462"""

    def __init__(self, pair_usage_table=None, location=None, species: str = 'h_sapiens'):
        super().__init__()

        if pair_usage_table is None:
            with open(os.path.join(BACKEND_OBJECTIVES_DATA, f'dinucleotides_usage_{species}.json'),
                      'r') as myfile:
                data = myfile.read()
            self.pair_usage_table = json.loads(data)
        else:
            self.pair_usage_table = pair_usage_table

        self.localization_group_spread = 5
        self.location = Location.from_data(location)

    def evaluate(self, problem: DnaOptimizationProblem) -> SpecEvaluation:
        """Return discrepancy score for the problem"""
        sequence = problem.sequence
        calculated_freqs = self.calculate_freqs(sequence)
        score = self.calculate_score(calculated_freqs)
        breaches = self.get_overrepresented_dinucleotides(sequence, calculated_freqs)

        return SpecEvaluation(
            self, problem,
            score=-score,
            locations=self.positions_to_locations(breaches),
            message=f"Dinucleotide score: {score}"
        )

    def __str__(self) -> str:
        """String representation."""
        return "MatchTargetPairUsage"

    def calculate_freqs(self, sequence: str) -> Dict[str, float]:
        """
        Calculate a dictionary of frequencies for each possible dinucleotide
        """
        total_number_of_dinuc = len(sequence) - 1
        calculated_freqs = {}

        counts = {}
        for i in range(0, len(sequence) - 1):
            counts[sequence[i:i + 2]] = counts.get(sequence[i:i + 2], 0) + 1
        for pair, _ in self.pair_usage_table.items():
            calculated_freqs[pair] = counts.get(pair, 0) / total_number_of_dinuc

        return calculated_freqs

    def calculate_score(self, calculated_freqs: Dict[str, float]) -> float:
        ratios = []
        for pair, fr in calculated_freqs.items():
            ratios.append(fr / self.pair_usage_table[pair])
        return mean([abs(1 - ratio) for ratio in ratios])

    def get_overrepresented_dinucleotides(self, sequence: str, freqs: Dict[str, float]) -> List[int]:
        dinuc_positions = {dinuc: [] for dinuc in self.pair_usage_table.keys()}
        for i in range(len(sequence) - 1):
            dinucleotide = sequence[i:i + 2]
            dinuc_positions[dinucleotide].append(i)
        nonoptimal = []
        for dinuc, freq in freqs.items():
            if freq > self.pair_usage_table[dinuc]:
                nonoptimal += dinuc_positions[dinuc]
        return nonoptimal

    def positions_to_locations(self, positions: List[int]) -> List[Location]:
        groups = group_nearby_indices(positions, max_group_spread=self.localization_group_spread)
        return [
            Location(group[0], group[-1] + 2, strand=1)
            for group in groups
        ]

    def localized(self, location, problem=None):
        return self.copy_with_changes(location=location)

```
    
# common/objectives/MFE.py
```python
from common.utils.Logger import MyLogger

logger = MyLogger(__name__)


class MFE():
    """This class implements the MFE prediction using a correlated stem–loop prediction.
    The method is described in Nucleic Acids Research, 2013, Vol. 41, No. 6 e73: 'mRNA secondary 
    structure optimization using a correlated stem–loop prediction'  """
    def __init__(self):
        pass

    def get_energy(self, seq1: str, seq2: str) -> float:
        """Calculation of contribution of each binding pair"""

        bond_energy = 0
        for i in range(len(seq1)):
            try:
                if [seq1[i], seq2[i]] in [['G', 'C'], ['C', 'G']]:
                    e = 3.12
                elif [seq1[i], seq2[i]] in [['A', 'T'], ['T', 'A']]:
                    e = 1
                elif [seq1[i], seq2[i]] in [['G', 'T'], ['T', 'G']]:
                    e = 1
                else:
                    e = 0
            except Exception as exception:
                #TODO: Explore why this exception happens 
                if str(exception) == "string index out of range":
                    e = 0
                    logger.error(f'{str(exception)}')
                    logger.debug(f'i: {i}, len(seq1): {len(seq1)}, len(seq2): {len(seq2)}')
                else:
                    #TODO: implement proper error handling
                    raise ValueError(str(exception))
            bond_energy = bond_energy + e

        return bond_energy

    def estimate(self, sequence: str) -> float:
        """Main method to get the estimation of the MFE"""

        sequence_size = len(sequence)
        iblock_size = 2 # <-- initial block size
        fblock_size = sequence_size/2 # <-- final block size
        loop_size = 3 # <-- minimum loop size
        c_energy = 0 # <-- cumulative energy

        for i in range(2):
            b = iblock_size
            while (b < fblock_size) and (sequence_size > loop_size + 2*b + 1): #  <-- Beware of never ending loop with celery!
                b += 1
                sub_seq1 = sequence[0:b]
                sub_seq2 = sequence[loop_size + b:loop_size + 2*b]
                energy = self.get_energy(sub_seq1, sub_seq2)
                c_energy = c_energy + energy
            sequence = sequence[::-1]

        return -c_energy/sequence_size

```
    
# common/objectives/MFE_optimization.py
```python
from common.objectives.MFE import MFE
from dnachisel import (Specification, SpecEvaluation, Location, DnaOptimizationProblem)
from common.utils.Logger import MyLogger
import RNA

logger = MyLogger(__name__)


class MinimizeMFE(Specification):
    """Minimize absolute value of MFE in the entropy window"""

    def __init__(self, five_end: str, entropy_window: int, boost: float = 10.0, mfe_method: str = "stem-loop") -> None:
        super().__init__()
        self.boost = boost
        self.five_end = five_end
        self.mfe_method = mfe_method
        self.entropy_window = entropy_window

    def localized(self, location, problem=None):
        return self.copy_with_changes(location=location)

    def evaluate(self, problem: DnaOptimizationProblem) -> SpecEvaluation:
        """Return MFE score for the problem' sequence in entropy window"""
        sequence = problem.sequence
        mfe_estimator = MFE()
        loc = Location(0, self.entropy_window)
        try:
            if self.mfe_method == "stem-loop":
                mfe = mfe_estimator.estimate(self.five_end + sequence[:self.entropy_window])
            elif self.mfe_method == "RNAfold":
                _, mfe = RNA.fold(self.five_end + sequence[:self.entropy_window])
            else:
                logger.warning(f'MFE method {self.mfe_method} is unknown! Switching to default one')
                _, mfe = RNA.fold(self.five_end + sequence[:self.entropy_window])
        except Exception as e:
            logger.error(f'Error when calling MFE.estimate {str(e)}')
        score = 100*mfe
        return SpecEvaluation(
            self, problem,
            score=score,
            locations=[loc],
            message="MFE score: {}".format(score)
        )

    def __str__(self) -> str:
        """String representation."""
        return "MinimizeMFE"

```
    
# common/objectives/Objectives.py
```python
from typing import List, Union, Tuple

from dnachisel import *
from common.objectives.Codon_pair_usage import MatchTargetCodonPairUsage
from common.objectives.Dinucleotide_usage import MatchTargetPairUsage
from common.objectives.MFE_optimization import MinimizeMFE
from common.utils.Logger import MyLogger
from common.utils.Exceptions import UnsupportedOptimizationError

# Setting up a logger
logger = MyLogger(__name__)


class Objectives(object):

    def __init__(self, entropy_window: int, organism: str, location: Union[None, Tuple], min_GC_content: int,
                 max_GC_content: int, GC_window_size: int, five_end: str, mfe_method: str, optimization_criterion: str) -> None:
        self.entropy_window = entropy_window
        self.organism = organism
        self.location = location
        self.min_GC_content = min_GC_content
        self.max_GC_content = max_GC_content
        self.GC_window_size = GC_window_size
        self.five_end = five_end
        self.mfe_method = mfe_method
        self.optimization_criterion = optimization_criterion

    def create_objectives(self) -> List[Specification]:
        """
        Create a list of objectives based on parameters.
        MinimizeMFE is added always, the others only if the corresponding parameter is set.
        MatchTargetCodonUsage is added when neither dinucleotides nor codon_pair is set. The three are mutually exclusive.
        """
        objectives = []
        objectives.append(MinimizeMFE(self.five_end, self.entropy_window, mfe_method=self.mfe_method))

        if self.GC_window_size:
            objectives.append(
                EnforceGCContent(mini=self.min_GC_content, maxi=self.max_GC_content, window=self.GC_window_size,
                                 location=self.location))


        if self.optimization_criterion == 'codon_usage':
            objectives.append(MatchTargetCodonUsage(species=self.organism, location=self.location))
        elif self.optimization_criterion == 'cai':
            objectives.append(MaximizeCAI(species=self.organism))
        elif self.optimization_criterion == 'dinucleotides':
            objectives.append(MatchTargetPairUsage(species=self.organism))
        elif self.optimization_criterion == 'codon_pair':
            objectives.append(MatchTargetCodonPairUsage(species=self.organism))
        else:
            logger.error(f'Unknown optimization criterion: {self.optimization_criterion}')
            raise UnsupportedOptimizationError(f'Unknown optimization criterion: {self.optimization_criterion}')

        return objectives

```
    
# common/objectives/__init__.py
```python

```
    
# common/utils/Datatypes.py
```python
from pydantic import BaseModel
from typing import List, Tuple


class OptimizationParameters(BaseModel):
    """Class for keeping all the input parameters for the optimization problem"""
    input_mRNA: str
    input_DNA: str
    five_end: str
    three_end: str
    avoid_motifs: List[str]
    max_GC_content: float
    min_GC_content: float
    GC_window_size: int
    usage_threshold: float
    uridine_depletion: bool
    organism: str
    entropy_window: int
    number_of_sequences: int
    filename: str
    mfe_method: str
    optimization_criterion: str
    location: Tuple[int, int, int]


class SequenceProperties(BaseModel):
    """Class for keeping track of sequence properties after evaluation"""
    seqID: str
    RNASeq: str
    DNASeq: str
    RNA_structure: str
    A_ratio: float
    TU_ratio: float
    G_ratio: float
    C_ratio: float
    AT_ratio: float
    GC_ratio: float
    MFE: float
    MFE_5: float
    score: float
    CAI: float


class EvaluationResults(BaseModel):
    """Class collecting the evaluation results together"""
    input: SequenceProperties
    optimized: List[SequenceProperties]

```
    
# common/utils/Exceptions.py
```python
class SequenceLengthError(Exception):
    """ This exception is raised when the input sequences have incorrect length"""

    def __init__(self, sequence_type):
        self.sequence_type = sequence_type
        self.message = "{} sequence with size not multiple of 3".format(self.sequence_type)
        super().__init__(self.message)


class OptimizationFailedError(Exception):
    """This exception is raised when the number of attempts is exceeded and no optimized
    sequences have been returned"""

    def __init__(self):
        self.message = 'Not possible to optimize with given parameters. Please try other ones.'
        super().__init__(self.message)


class EmptySequenceError(Exception):
    """This exception is raised when the input RNA or five and three flanking sequences are empty"""

    def __init__(self, sequence_type):
        self.sequence_type = sequence_type
        self.message = '{} sequence is not provided!'.format(sequence_type)
        super().__init__(self.message)


class WrongCharSequenceError(Exception):
    """This exception is raised when the input RNA or five and three flanking sequences contains forbidden characters"""

    def __init__(self, sequence_type):
        self.message = '{} sequence can only have "A", "T", "U", "G", "C" characters!'.format(sequence_type)
        super().__init__(self.message)


class RangeError(Exception):
    """This exception is raised when the min, max GC content or codon usage frequency threshold exceeds allowed range"""

    def __init__(self, message):
        super().__init__(message)


class EntropyWindowError(Exception):
    """This exception is raised when the entropy window is not defined by user"""

    def __init__(self):
        self.message = 'Entropy window is negative or not provided!'
        super().__init__(self.message)


class NoGCError(Exception):
    """This error is raised when the max or min GC is not provided by user"""

    def __init__(self):
        self.message = 'Max or Min GC content is not provided!'
        super().__init__(self.message)


class SpeciesError(Exception):
    """This error is raised when an invalid species is provided by user"""

    def __init__(self, message):
        super().__init__(message)


class NumberOfSequencesError(Exception):
    """This error is raised when user didn't input the number of sequences they
    want to see in the output """

    def __init__(self):
        self.message = 'Number of sequences in the output is not defined or negative!'
        super().__init__(self.message)


class UnsupportedOptimizationError(Exception):
    """This error is raised when unsupported criterion for codon optimization"""

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

```
    
# common/utils/Logger.py
```python
from loguru import logger
import os


def MyLogger(name):
    """ Logger which can be imported to any module """

    logger.remove()

    logger.level("DEBUG", color="<yellow>")
    logger.level("INFO", color="<green>")
    logger.level("ERROR", color="<red>")

    file = os.environ.get('LOG_FILE')

    if file:
        log_format = "{time:YYYY-MMMM-DD,dddd    HH:mm:ss:A}| {module}| {level}| {message}"
        logger.add(file, format=log_format, level="DEBUG", backtrace=True, rotation="100 MB", compression="zip",
                   enqueue=True)
    logger.enable(name)

    return logger

```
    
# common/utils/RequestParser.py
```python
import json

from common.utils.Datatypes import OptimizationParameters
from common.utils.Exceptions import SequenceLengthError, EmptySequenceError, \
    EntropyWindowError, NoGCError, NumberOfSequencesError, WrongCharSequenceError, RangeError, SpeciesError
from common.utils.Logger import MyLogger

# Setting up a logger
logger = MyLogger(__name__)


class RequestParser(object):
    """This class parses the input request and returns a set of
    parameters necessary to perform the optimization task """

    def __init__(self, request):
        self.data = json.loads(request.get_data())
        logger.info(self.data)

    def parse(self) -> OptimizationParameters:
        print(self.data)
        def seq_char_check(seq, seq_name):
            """Function to check if allowed characters are used"""
            if not set(seq) <= set("TAGC"):
                raise WrongCharSequenceError(seq_name)

        # --- Get the input sequences ---
        input_mRNA = self.data['sequences']['gene_of_interest'].upper().replace('T', 'U')
        input_DNA = input_mRNA.replace('U', 'T')
        if not input_mRNA:
            raise EmptySequenceError('RNA')
        if len(input_DNA) % 3:
            raise SequenceLengthError('Input RNA')
        seq_char_check(input_DNA, seq_name="RNA")
        five_end = self.data['sequences']['five_end_flanking_sequence'].upper().replace('U', 'T')
        if not five_end:
            raise EmptySequenceError('Five prime flanking')
        seq_char_check(five_end, seq_name="Five prime flanking")
        three_end = self.data['sequences']['three_end_flanking_sequence'].upper().replace('U', 'T')
        if not three_end:
            raise EmptySequenceError('Three prime flanking')
        seq_char_check(three_end, seq_name="Three prime flanking")

        # --- Get the list of avoided codons ---
        avoid_motifs = self.data['config']['avoided_motifs']

        # --- Get GC related parameters ---
        max_GC_content = self.data['config']['max_GC_content']
        min_GC_content = self.data['config']['min_GC_content']
        GC_window_size = self.data['config']['GC_window_size']
        if max_GC_content is None:
            logger.error('No max GC content defined')
            raise NoGCError()
        elif min_GC_content is None:
            logger.error('No min GC content defined')
            raise NoGCError()
        elif not GC_window_size:
            logger.warning('No GC window size defined or set to 0. Window GC optimization will not be performed.')
        elif max_GC_content > 1 or max_GC_content < 0:
            raise RangeError('Max GC content can be only in a range of 0 to 1!')
        elif min_GC_content > 1 or min_GC_content < 0:
            raise RangeError('Min GC content content can be only in a range of 0 to 1!')
        elif min_GC_content > max_GC_content:
            raise RangeError("Min GC content connot be higher than max GC content")
        else:
            pass

        # --- Get the threshold for usage frequency of codons ---
        usage_threshold = self.data['config']['codon_usage_frequency_threshold']
        if usage_threshold is None:
            logger.warning('Usage threshold is not defined')
        elif usage_threshold > 1 or usage_threshold < 0:
            raise RangeError('Codon usage frequency threshold can be only in a range of 0 to 1!')
        else:
            pass

        # --- Get the uridine depletion boolean value ---
        uridine_depletion = self.data['uridine_depletion']

        # --- Get the organism ---
        organism = self.data['config']['organism']
        if organism.lower() in ["h_sapiens", "sapiens", "human", "homo_sapiens", "h sapiens", "homo sapiens"]:
            organism = "h_sapiens"
        elif organism.lower() in ["m_musculus", "mouse", "mus_musc", "mus_musculus", "mus", "mus musculus", "m musculus", "mus musc"]:
            organism = "m_musculus"
        elif organism.lower() in ["custom"]:
            organism = "h_sapiens"
        else:
            raise SpeciesError("Invalid organism {} is provided!".format(organism))

        # --- Get the entropy window value ---
        entropy_window = self.data['config']['entropy_window']
        if not entropy_window:
            logger.warning('Entropy window is not defined or is zero')
            raise EntropyWindowError()
        elif entropy_window < 0:
            logger.warning('Entropy window is negative, meaning is ambiguous.')
            raise EntropyWindowError()

        # --- Get the number of sequences required by user ---
        number_of_sequences = self.data['config']['number_of_sequences']
        if not number_of_sequences or number_of_sequences < 0:
            logger.error('Number of sequences for the output is not defined or negative!')
            raise NumberOfSequencesError()

        # --- Get the filename to produce the report ---
        filename = self.data['file_name']
        if not filename:
            filename = 'optimization_results'

        # --- Get the algorithm for MFE estimation ---
        precise_MFE_algorithm = self.data['precise_MFE_algorithm']
        if precise_MFE_algorithm is None:
            logger.warning('Algorithm for MFE is not specified! Switching to the default one...')
            precise_MFE_algorithm = False

        if precise_MFE_algorithm:
            mfe_method = 'RNAfold'
        else:
            mfe_method = 'stem-loop'

        # --- Get the optimization criterion string value
        optimization_criterion = self.data['optimization_criterion']

        # --- Combining all the parameters together ---
        return OptimizationParameters(input_mRNA=input_mRNA,
                                      input_DNA=input_DNA,
                                      five_end=five_end,
                                      three_end=three_end,
                                      avoid_motifs=avoid_motifs,
                                      max_GC_content=max_GC_content,
                                      min_GC_content=min_GC_content,
                                      GC_window_size=GC_window_size,
                                      usage_threshold=usage_threshold,
                                      uridine_depletion=uridine_depletion,
                                      organism=organism,
                                      entropy_window=entropy_window,
                                      number_of_sequences=number_of_sequences,
                                      filename=filename,
                                      mfe_method=mfe_method,
                                      optimization_criterion=optimization_criterion,
                                      location=(0, len(input_mRNA) - len(input_mRNA) % 3, 1))

```
    
# common/utils/__init__.py
```python

```
    
# tests/test_CAI.py
```python
from common.Evaluation import Evaluation
from common.utils.Datatypes import OptimizationParameters
from unittest import main, TestCase

test_cases = [('ATGAAGTCGCAATTCAGCTCGCACGTACTCCCGCTAAGCCGAAAGGTCCTATAATCTAGCCGGTCGGGTGTCTGAGATAACTATGCTAACCAGGCCTAACCGCGCGGTATTCAGCTAAGCACTCCGATCCCAAGAGCCAACCGCGAACCAAGTGAAGACTATATCTGTTCGTTCTTCACG', 0.63),
              ('GTAATAGGACTTCTCAAGCAGGTTATGTACCGTCATGTCTTGACCTCCTACGCTTAATCTCTTAGGGGCCCACGGACCTGCGTTCTTCGTCCGCATGACTTACTGCCACTCATCTGTCCGTCTCCTCCCCAGTTGAACCCACGATCGTGG', 0.62),
              ('CCGAAACAATTGGGACACTAGATACAGGTTAACTCCATGGCAGAGACCAATAAACGACCCTGGAAGATTGGCAGCCTTTGAAAATACACGGCGCGTTATACCTCTCAACCTCTGCTGCGCCATTGGCCAGATAGGTCCCCTCCTCTCAGCCCTAGGGCAAAGCTCGAACGTACTCTAGGTTACCCTTACTGCACTTAGATATTCTTCTCCCGGTGACTTGCCGATGATTCTACAAGGGTCGGCTGGGTCGATCCACGTGGTTACCGTGAT', 0.7),
              ('CTAAATGTGATTCCAGTGCGTTGATCGTCGGTTCCATTCACTCTCTCGCTCTATAGTGAGGTTGGTCATCGCATATAACTGCTACTCTGCCTATTGTACAAAACTATCTTGTTATGTTGCGACAACGCCAGATAGCATGTTGGTTTATAGACGTCTTGTTTGCGGGCTAAGGGATACAGAACTGTTGAATGATATCAATCTGCGATCGAGTCGCTTCATAAAGTGTGCGTCTGTCGACTGAGCTGTGGTCGCTGGATGCTCAGTTCTGGT', 0.6)]

class TestCAI(TestCase):

    def test_h_sapiens(self):
        parameters = OptimizationParameters(
            organism='h_sapiens',
        )
        evaluation = Evaluation(parameters=parameters)
        for seq in test_cases:
            self.assertEqual(round(evaluation.CAI(seq[0]), 2), seq[1])
            
if __name__ == '__main__':
    main()
```
    
# tests/test_codon_pair.py
```python
from statistics import mean
from unittest import TestCase, main
from dnachisel import DnaOptimizationProblem
from common.objectives.Codon_pair_usage import MatchTargetCodonPairUsage

def get_ratios_from_sequence(sequence):
    """Function to calculate codon pair frequencies from the input sequence"""
    objective = MatchTargetCodonPairUsage()
    pair_usage_table = objective.pair_usage_table
    total_number_of_pairs = len(sequence) / 3 - 1
    calculated_freqs = {}
    counts = {}
    for i in range(0, len(sequence) - 3, 3):
        counts[sequence[i:i + 6]] = counts.get(sequence[i:i + 6], 0) + 1
    for pair, _ in counts.items():
        calculated_freqs[pair] = counts.get(pair, 0) / total_number_of_pairs
    return calculated_freqs


def calculate_mean_ratio(calculated_freqs):
    """Function which calculates mean value from the dictionary of frequencies"""
    objective = MatchTargetCodonPairUsage()
    pair_usage_table = objective.pair_usage_table
    ratios = []
    for pair, fr in calculated_freqs.items():
        ratios.append(pair_usage_table[pair] / fr)
    return mean([abs(1 - ratio) for ratio in ratios])

class TestCodonPair(TestCase):
    def optimization(self, organism):
        """Check that optimized codon pair frequencies are closer to target than original ones"""
        input_sequence = "CCGTCGCGGCAGGTTATTATACCTCATTCCTTGGAGACATACAACTATCAATGGGACTTGAGGTTAAGGTATTCCCGCATGAACGCGTGTACTGAAAATATGAAGGCGAGGGCGGAAGCTTTCATTAGCGAGCACCTACAACGTTAGAGTTGGTCGTGTCTTGCTATGCGTCCAGCACATCTGTAAGCCGGTATAAGGCCAGGGGCGGTACATATCGTACAGATCTAGTACATGTTGATAACTTTCATCTGTCGTAGGAAGGCGGAGCCGCCCCTGACGGACGTAGAAAGGGGAATGGGCACTGAGACCCAGTGAGCCCCTTTTGCGTTCTTGGCAAATACCTAGACCTTCTGGTCGTCCTATCGTAATATCTCCTGATACTCATGACAGCAGGATAGCAGCCTGCAACCTCCATGTACTTCGTTGGATTCTTTCCGAGTCTCGTGTGAGTAGATGCTTTGGGGAGTTACCTCTAACACATGGCTTGTTTATTCGTAATTCGACTCCCATGCTTGCTTTTAAACGTCTGTCAACATGAACATTCTGGTCGCACGACGATTAAGAAAGGGAACTTCGTGTTGATGTAGTAGGATATAGCAG"
        objective = MatchTargetCodonPairUsage(species=organism)
        problem = DnaOptimizationProblem(input_sequence, objectives=[objective])
        problem.resolve_constraints()
        problem.optimize()

        optimized_sequence = problem.sequence
        optimized_ratios = get_ratios_from_sequence(optimized_sequence)
        optimized_mean = calculate_mean_ratio(optimized_ratios)

        input_ratios = get_ratios_from_sequence(input_sequence)
        input_mean = calculate_mean_ratio(input_ratios)

        assert optimized_mean < input_mean

    def test_optimized_frequencies_human(self):
        self.optimization('h_sapiens')

    def test_optimized_frequencies_mouse(self):
        self.optimization('m_musculus')
```
    
# tests/test_dinucleotide_frequencies.py
```python
from statistics import mean

from Bio.Seq import Seq
from dnachisel import DnaOptimizationProblem, EnforceTranslation, Location
from common.objectives.Dinucleotide_usage import MatchTargetPairUsage


def test_frequencies():
    """Check calculation of frequencies"""
    seq1 = "AGCAGAGAAGGCGGAAGCAGTGGCGTCCGCAGCTGGGGCTT"
    """
    AG GC CA AG GA AG GA AA AG GG GC CG GG GA AA AG GC CA AG GT TG GG GC CG GT TC CC CG GC CA AG GC CT TG GG GG GG GC CT TT
    40 dinucleotides  
    """
    counts = {"TT": 1, "TC": 1, "TA": 0, "TG": 2,
              "CT": 2, "CC": 1, "CA": 3, "CG": 3,
              "AT": 0, "AC": 0, "AA": 2, "AG": 7,
              "GT": 2, "GC": 7, "GA": 3, "GG": 6}

    total_count = 40
    objective = MatchTargetPairUsage()

    calculated_freqs = objective.calculate_freqs(seq1)
    for pair, value in calculated_freqs.items():
        assert value == counts[pair] / total_count
    score = objective.calculate_score(calculated_freqs)
    print(score)


def get_ratios_from_sequence(sequence):
    """Function to calculate dinucleotide frequencies from the input sequence"""
    objective = MatchTargetPairUsage()
    pair_usage_table = objective.pair_usage_table
    total_number_of_dinuc = len(sequence) - 1
    calculated_freqs = {}
    counts = {}
    for i in range(0, len(sequence) - 1):
        counts[sequence[i:i + 2]] = counts.get(sequence[i:i + 2], 0) + 1
    for pair, _ in pair_usage_table.items():
        calculated_freqs[pair] = counts.get(pair, 0) / total_number_of_dinuc
    return calculated_freqs


def calculate_mean_ratio(calculated_freqs):
    """Function which calculates mean value from the dictionary of frequencies"""
    objective = MatchTargetPairUsage()
    pair_usage_table = objective.pair_usage_table
    ratios = []
    for pair, fr in calculated_freqs.items():
        ratios.append(fr / pair_usage_table[pair])
    return mean([abs(1 - ratio) for ratio in ratios])

def optimization(organism: str):
    """Check that optimized dinucleotide frequencies are closer to target than original ones"""
    input_sequence = "CCGTCGCGGCAGGTTATTATACCTCATTCCTTGGAGACATACAACTATCAATGGGACTTGAGGTTAAGGTATTCCCGCATGAACGCGTGTACTGAAAATATGAAGGCGAGGGCGGAAGCTTTCATTAGCGAGCACCTACAACGTTAGAGTTGGTCGTGTCTTGCTATGCGTCCAGCACATCTGTAAGCCGGTATAAGGCCAGGGGCGGTACATATCGTACAGATCTAGTACATGTTGATAACTTTCATCTGTCGTAGGAAGGCGGAGCCGCCCCTGACGGACGTAGAAAGGGGAATGGGCACTGAGACCCAGTGAGCCCCTTTTGCGTTCTTGGCAAATACCTAGACCTTCTGGTCGTCCTATCGTAATATCTCCTGATACTCATGACAGCAGGATAGCAGCCTGCAACCTCCATGTACTTCGTTGGATTCTTTCCGAGTCTCGTGTGAGTAGATGCTTTGGGGAGTTACCTCTAACACATGGCTTGTTTATTCGTAATTCGACTCCCATGCTTGCTTTTAAACGTCTGTCAACATGAACATTCTGGTCGCACGACGATTAAGAAAGGGAACTTCGTGTTGATGTAGTAGGATATAGCAG"
    objective = MatchTargetPairUsage(species=organism)
    problem = DnaOptimizationProblem(input_sequence, objectives=[objective])
    problem.resolve_constraints()
    problem.optimize()

    optimized_sequence = problem.sequence
    # optimized_sequence = "CCAGCCCGGGTCGCCGGCTACGCAGGTTATTGTTGGACACAAGCAACAAACTAGATCGGAAGGTCTTATATCGTAATACAGACTAGGAATTTGGCGCGCTGGCGTCTGCTTCAAGTCCGCCAGGGTAATGGCCTTGAATACGGATCCCCAGACCTTCTCGCGACCATGCCAGGGAGAGAGCAGTCCATGCTAAACTCCGACCCACGGTCAACATCAATCCGTGTAAGTAGGGTTATTCATAACTATATCAGCAGCGTGGCCGATGGCGACAATTTAGGTCTTAAGTCGGGAGGCAACACCTTATACCCTACATCGTAGACCGCCGCCAACGACGTAATTTGGCTGTCAAGGAGGCGGGCGTGTGAGACCTCTTGTTCGTCTATAGAGACACGGTCATCCCATATTTTATGTAGATCGCGACGTGAAATTGAGCTAATTGGCCGTCTTGATTGAAAGGGAGAACCCTATCCTGTGGGAATTACAAATCGGTCTTGATTGTGGAAGGTTCCTAGGAATCATGCCTCAAGCGCATACAGTAGTGCCGCACATTAAACCGTAGTGTCTAGGATGCTTACGGGGCTGACTGGCCTCATGGGATCT"
    optimized_ratios = get_ratios_from_sequence(optimized_sequence)
    optimized_mean = calculate_mean_ratio(optimized_ratios)

    input_ratios = get_ratios_from_sequence(input_sequence)
    input_mean = calculate_mean_ratio(input_ratios)

    assert optimized_mean < input_mean


def test_optimized_frequencies_human():
    optimization('h_sapiens')

def test_optimized_frequencies_mouse():
    optimization('m_musculus')



def test_correct_translation():
    input_sequence = "AUGACGAAAUCAGAUAUAGCUAAGGAACUCGCAAGGAGGCACGGCAUAUCCUACAAAAAAGCCCUCCUUAUAGUCAACAUGACCUUUGAGAUACUAAAAGCAAAAAUACUGAACGGUGAAAAGGUAGAGGUAAGGGGACUGGGAACCUUUAAGUUGAAGAGAAAACCGGGAAGGUUCGUAAAGAACCCGAAAACGGGUAUAGAAAUUUACGUAAAGGAGAGGUACGUUCCCUACUACAAGAUGUCUAAACUUUUAAGGAAGAAACUAAACGGCGAUAAAGAAAGGGAGGAGUGUUUGACUUGA"
    original_translation = Seq(input_sequence).translate()
    objective = MatchTargetPairUsage(species='h_sapiens')
    constraint = EnforceTranslation(genetic_table='Standard', location=Location(0, len(input_sequence), 1))

    problem = DnaOptimizationProblem(input_sequence.replace('U', 'T'), objectives=[objective], constraints=[constraint])
    problem.resolve_constraints()
    problem.optimize()
    optimized_sequence = problem.sequence
    translation = Seq(optimized_sequence).translate()

    assert translation == original_translation

    optimized_ratios = get_ratios_from_sequence(optimized_sequence)
    optimized_mean = calculate_mean_ratio(optimized_ratios)

    input_ratios = get_ratios_from_sequence(input_sequence)
    input_mean = calculate_mean_ratio(input_ratios)

    assert optimized_mean < input_mean

```
    
# tests/test_exception_empty_seq.py
```python
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

```
    
# tests/test_exception_entropy_win.py
```python
from tests.test_request_parsing import MockRequest
from common.utils.Exceptions import EntropyWindowError
from common.utils.RequestParser import RequestParser


def test_neg_entropy_window():
    config = {
        "avoid_codons": ["ATG"],
        "avoided_motifs": [],
        "codon_usage_frequency_threshold": 0.1,
        "max_GC_content": 0.7,
        "GC_window_size": 30,
        "min_GC_content": 0.3,
        "organism": 'h_sapiens',
        "entropy_window": -6,
        "number_of_sequences": 5
    }
    dinucleotides = False
    uridine_depletion = True
    faster_MFE_algorithm = True
    filename = 'test'
    sequences = {"five_end_flanking_sequence": "ACG",
                 "gene_of_interest": "AGCAAAGAAGGCGGA",
                 "three_end_flanking_sequence": "ACG"}
    request = MockRequest(config=config,
                          dinucleotides=dinucleotides, uridine_depletion=uridine_depletion,
                          faster_MFE_algorithm=faster_MFE_algorithm, filename=filename, sequences=sequences)
    parser = RequestParser(request)

    try:
        parameters = parser.parse()
        right_error_raised = False
    except EntropyWindowError as e:
        right_error_raised = True
    except Exception as e:
        right_error_raised = False
    assert right_error_raised


def test_no_entropy_window():
    config = {
        "avoid_codons": ["ATG"],
        "avoided_motifs": [],
        "codon_usage_frequency_threshold": 0.1,
        "max_GC_content": 0.7,
        "GC_window_size": 30,
        "min_GC_content": 0.3,
        "organism": 'h_sapiens',
        "entropy_window": None,
        "number_of_sequences": 5
    }
    dinucleotides = False
    uridine_depletion = True
    faster_MFE_algorithm = True
    filename = 'test'
    sequences = {"five_end_flanking_sequence": "ACG",
                 "gene_of_interest": "AGCAAAGAAGGCGGA",
                 "three_end_flanking_sequence": "ACG"}
    request = MockRequest(config=config,
                          dinucleotides=dinucleotides, uridine_depletion=uridine_depletion,
                          faster_MFE_algorithm=faster_MFE_algorithm, filename=filename, sequences=sequences)
    parser = RequestParser(request)

    try:
        parameters = parser.parse()
        right_error_raised = False
    except EntropyWindowError as e:
        right_error_raised = True
    except Exception as e:
        right_error_raised = False
    assert right_error_raised

```
    
# tests/test_exception_missing_parameter.py
```python
from tests.test_request_parsing import MockRequest
from common.utils.RequestParser import RequestParser


def test_parse_incomplete_parameters():
    config = {
        "avoid_codons": [],
        "max_GC_content": 0.7,
        "GC_window_size": 30,
        "min_GC_content": 0.3,
        "entropy_window": 6,
        "number_of_sequences": 5
    }
    sequences = {
        "five_end_flanking_sequence": "AAA",
        "gene_of_interest": "AAAGGGCCCTTT",
        "three_end_flanking_sequence": "AAA"
    }
    request = MockRequest(config=config, sequences=sequences, faster_MFE_algorithm=False)

    parser = RequestParser(request)
    try:
        parameters = parser.parse()
        right_error_raised = False
    except KeyError as e:
        right_error_raised = True
    except:
        right_error_raised = False
    assert right_error_raised

```
    
# tests/test_exception_no_gc.py
```python
from tests.test_request_parsing import MockRequest
from common.utils.Exceptions import NoGCError
from common.utils.RequestParser import RequestParser


def test_no_gc():
    config = {
        "avoid_codons": ["ATG"],
        "avoided_motifs": [],
        "codon_usage_frequency_threshold": 0.1,
        "max_GC_content": None,
        "GC_window_size": 30,
        "min_GC_content": 0.3,
        "organism": 'h_sapiens',
        "entropy_window": 10,
        "number_of_sequences": 5
    }
    dinucleotides = False
    uridine_depletion = True
    faster_MFE_algorithm = True
    filename = 'test'
    sequences = {"five_end_flanking_sequence": "ACG",
                 "gene_of_interest": "AGCAAAGAAGGCGGA",
                 "three_end_flanking_sequence": "ACG"}
    request = MockRequest(config=config,
                          dinucleotides=dinucleotides, uridine_depletion=uridine_depletion,
                          faster_MFE_algorithm=faster_MFE_algorithm, filename=filename, sequences=sequences)
    parser = RequestParser(request)
    try:
        parameters = parser.parse()
    except NoGCError as e:
        right_error_raised = True
    except Exception as e:
        print(e)
        right_error_raised = False
    else:
        right_error_raised = False
    assert right_error_raised

```
    
# tests/test_exception_num_of_seqs.py
```python
from tests.test_request_parsing import MockRequest
from common.utils.Exceptions import NumberOfSequencesError
from common.utils.RequestParser import RequestParser


def test_num_of_seqs():
    config = {
        "avoid_codons": ["ATG"],
        "avoided_motifs": [],
        "codon_usage_frequency_threshold": 0.1,
        "max_GC_content": 0.7,
        "GC_window_size": 30,
        "min_GC_content": 0.3,
        "organism": 'h_sapiens',
        "entropy_window": 10,
        "number_of_sequences": 0
    }
    dinucleotides = False
    uridine_depletion = True
    faster_MFE_algorithm = True
    filename = 'test'
    sequences = {"five_end_flanking_sequence": "ACG",
                 "gene_of_interest": "AGCAAAGAAGGCGGA",
                 "three_end_flanking_sequence": "ACG"}
    request = MockRequest(config=config,
                          dinucleotides=dinucleotides, uridine_depletion=uridine_depletion,
                          faster_MFE_algorithm=faster_MFE_algorithm, filename=filename, sequences=sequences)
    parser = RequestParser(request)

    try:
        parameters = parser.parse()
        right_error_raised = False
    except NumberOfSequencesError as e:
        right_error_raised = True
    except KeyError as e:
        right_error_raised = True
    assert right_error_raised

```
    
# tests/test_exception_optimization_failure.py
```python
import pytest
from Bio.Data.CodonTable import TranslationError
from common.utils.Exceptions import OptimizationFailedError
from common.OptimizationProblems import initialize_optimization_problem
from common.OptimizationTask import optimization_task
from common.utils.Datatypes import OptimizationParameters


@pytest.fixture
def example_parameters():
    def _example_parameters(mRNA_seq):
        dna = mRNA_seq.replace('U', 'T')
        parameters = OptimizationParameters(
            input_mRNA=mRNA_seq,
            input_DNA=dna,
            five_end="ATG",
            three_end="CGG",
            avoid_codons=['AAA', 'AAG'],
            avoid_motifs=[],
            max_GC_content=0.9,
            min_GC_content=0.8,
            GC_window_size=30,
            usage_threshold=0.1,
            uridine_depletion=True,
            organism='h_sapiens',
            entropy_window=8,
            number_of_sequences=3,
            mfe_method='stem_loop',
            dinucleotides=False,
            codon_pair=False,
            CAI=False,
            location=(0, len(mRNA_seq) - len(mRNA_seq) % 3, 1),
            filename="optimization_results"
        )
        return parameters
    return _example_parameters


def test_optimization_failure_1(example_parameters):
    """Assert that optimization fails when conflicting parameters are given"""
    seq = "AGCAAAGAAGGCGGA"
    params = example_parameters(mRNA_seq=seq)
    optimization_problem = initialize_optimization_problem(params)
    try:
        sequence = optimization_task(optimization_problem)
        right_error_raised = False
    except OptimizationFailedError as e:
        right_error_raised = True
    except:
        right_error_raised = False
    assert right_error_raised

```
    
# tests/test_exception_range_error.py
```python
from common.utils.Exceptions import RangeError
from common.utils.RequestParser import RequestParser
import pytest
from tests.test_request_parsing import MockRequest

def assert_range_error(request):
    parser = RequestParser(request)
    try:
        parameters = parser.parse()
    except RangeError as e:
        right_error_raised = True
    else:
        right_error_raised = False
    assert right_error_raised

@pytest.fixture
def example_parameters():
    def _example_parameters(max_GC_content, min_GC_content, codon_usage):
        return {
        "config": {
            "avoid_codons": [],
            "avoided_motifs": [],
            "codon_usage_frequency_threshold": codon_usage,
            "max_GC_content": max_GC_content,
            "GC_window_size": 30,
            "min_GC_content": min_GC_content,
            "organism": 'h_sapiens',
            "entropy_window": 6,
            "number_of_sequences": 5
        },
        "dinucleotides": False,
        "match_codon_pair": False,
        "uridine_depletion": True,
        "faster_MFE_algorithm": False,
        "precise_MFE_algorithm": False,
        "CAI": False,
        "file_name": "test",
        "sequences": {
            "five_end_flanking_sequence": "AGAACCGCACTCGGGAGGTTAAGCCGTCAATAACAATGACTCAAAGCCAAATGGTACTTGATAAGCGCAGCGTGGTGT",
            "gene_of_interest": "ATACTTTTCTGCATCGACTAAGTCTGAACGGAGAACCCGGAAACATATGTTGCTTCAGAACGCCGAGATTACTGCGTTTAGACTCAAACGAGAAGCGCCCCGCTCGTGTTATCTTCCAACTACCACGTATGTGTTGGACTGTTAGTGCACTTCAGATTGGTCCAAAACACCCTAGATGCTGAACCTGTACCCAGTAAGAAGTACGTAGTCGTTACAAGCAAACCCTAGCTTAGTAAAGGAAGCGCGTCTTACGTGTATGGTGGTTATCTTCCACCGTTGGACTTACGCGATTTGTGGAAATTACTTGCTTGCGTTGCTACGGGATTGAGAAAGCTTGCATATTGGGGGGTCAGACCCGTCATACAGGCCATCATGTCGGAAGGGCAGGTTCCATTGCCAGGTAGCTCCACCGCCCTGACGGGTTTCCAGGCACAGGTATACTCTCATACCAACAAAGCGCTTTCTGGGACTCGACGACACCAGGCCACGCGATGCTTCTCTAGATTCTCTGCGTTCGCTGTTGCATGCAGCCTCGCCTCGTTGAAGATGAGTATACAAACCCTGAACTGGGTGGTCCGATTGCTTGTAGGGTACACAACATGAGGGCGCACAAGCCCGGCGGGCCAATCCTTGATCACAACATATTCTACCTTCACTCCGCCAGACTGCGACCGCGTCAGTGTAAGCCCAAGTTTCGTGTTAGAAATACCACATGAATCGCGATTAACGGAGCACACGGCGCACGACGTCGATCCCTACTTAAGGGGTGTTTGACTCTCACAGTATGTTTATGCTCCGTACAGGCGCCCTAGTCGCAAGAAAGCCACTATGCTTGCTAATTTTTTCTTTGGCACGAACGTGCTTTCTAACTGGTAGATCCTATACTAGGCCTACTCGCTA",
            "three_end_flanking_sequence": "CATCGCGGGAATGCTACCTAGACTTGGGGTAAGTCGATCCAGTGTCCC"
        }
    }
    return _example_parameters


def test_max_GC_higher_1(example_parameters):
    """Assert that max GC content doen`t exceed 100%"""
    request = MockRequest(**example_parameters(1.2, 0.4, 0.1))
    assert_range_error(request)


def test_max_GC_less_0(example_parameters):
    """Assert that max GC content isn`t negative"""
    request = MockRequest(**example_parameters(-2.3, 0.4, 0.1))
    assert_range_error(request)


def test_min_GC_higher_1(example_parameters):
    """Assert that min GC content doen`t exceed 100%"""
    request = MockRequest(**example_parameters(0.7, 1.4, 0.1))
    assert_range_error(request)


def test_min_GC_less_0(example_parameters):
    """Assert that min GC content isn`t negative"""
    request = MockRequest(**example_parameters(0.7, -0.4, 0.1))
    assert_range_error(request)


def test_codon_usage_higher_1(example_parameters):
    """Assert that codon usage frequency threshold doen`t exceed 100%"""
    request = MockRequest(**example_parameters(0.7, 0.4, 3))
    assert_range_error(request)

def test_codon_usage_less_0(example_parameters):
    """Assert that codon usage frequency threshold isn`t negative"""
    request = MockRequest(**example_parameters(0.7, 0.4, -0.1))
    assert_range_error(request)


def test_max_GC_less_min_GC(example_parameters):
    """Assert that max GC content is higher than min GC content"""
    request = MockRequest(**example_parameters(0.1, 0.4, 0.1))
    assert_range_error(request)
```
    
# tests/test_exception_seq_length.py
```python
from tests.test_request_parsing import MockRequest
from common.utils.Exceptions import SequenceLengthError
from common.utils.RequestParser import RequestParser


def test_wrong_sequence_length():
    """Assert that SequenceLengthError is raised when sequence length is not divisible by 3"""
    config = {
        "avoid_codons": ["ATG"],
        "avoided_motifs": [],
        "codon_usage_frequency_threshold": 0.1,
        "max_GC_content": 0.7,
        "GC_window_size": 30,
        "min_GC_content": 0.3,
        "organism": 'h_sapiens',
        "entropy_window": 5,
        "number_of_sequences": 5
    }
    dinucleotides = False
    uridine_depletion = True
    faster_MFE_algorithm = True
    filename = 'test'
    sequences = {"five_end_flanking_sequence": "ACG",
                 "gene_of_interest": "AGCAAAGAAGGCGG",
                 "three_end_flanking_sequence": "ACG"}
    request = MockRequest(config=config,
                          dinucleotides=dinucleotides, uridine_depletion=uridine_depletion,
                          faster_MFE_algorithm=faster_MFE_algorithm, filename=filename, sequences=sequences)
    parser = RequestParser(request)
    try:
        parameters = parser.parse()
    except SequenceLengthError as e:
        right_error_raised = True
    except:
        right_error_raised = False
    else:
        right_error_raised = False
    assert right_error_raised

```
    
# tests/test_exception_species_error.py
```python
from common.OptimizationProblems import initialize_optimization_problem
from common.OptimizationTask import optimization_task
from common.utils.Datatypes import OptimizationParameters
from common.utils.Exceptions import SpeciesError


def test_species_error():
    mRNA_seq = "AAUCAAAUAGGGUUAAGUCUAGGAUUGUUAGUCUGCUAAGGUCUGCAGUUACUGUGUCUACUGAUGAUAGUUCGCAUUGACAAU"
    dna = mRNA_seq.replace('U', 'T')
    parameters = OptimizationParameters(
        input_mRNA=mRNA_seq,
        input_DNA=dna,
        five_end="ATG",
        three_end="CGG",
        avoid_codons=['AAA', 'AAG'],
        avoid_motifs=[],
        max_GC_content=0.7,
        min_GC_content=0.3,
        GC_window_size=30,
        usage_threshold=0.1,
        uridine_depletion=True,
        organism='wrong_organism',
        entropy_window=8,
        number_of_sequences=3,
        mfe_method='stem_loop',
        dinucleotides=False,
        CAI=False,
        codon_pair=False,
        location=(0, len(mRNA_seq) - len(mRNA_seq) % 3, 1),
        filename="optimization_results"
    )

    try:
        optimization_problem = initialize_optimization_problem(parameters)
        sequence = optimization_task(optimization_problem)
    except SpeciesError as e:
        right_error_raised = True
    except Exception as e:
        print(e)
        right_error_raised = False
    else:
        right_error_raised = False
    assert right_error_raised

```
    
# tests/test_exception_wrong_char_seq.py
```python
from common.utils.Exceptions import WrongCharSequenceError
from common.utils.RequestParser import RequestParser
import pytest
from tests.test_request_parsing import MockRequest

def assert_wrong_char(example_parameters, sequences):
    request = MockRequest(**example_parameters, sequences=sequences)
    parser = RequestParser(request)
    try:
        parameters = parser.parse()
    except WrongCharSequenceError as e:
        right_error_raised = True
    else:
        right_error_raised = False
    assert right_error_raised

@pytest.fixture
def example_parameters():
    return {
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
        "match_codon_pair": False,
        "uridine_depletion": True,
        "faster_MFE_algorithm": False,
        "precise_MFE_algorithm": False,
        "CAI": False,
        "file_name": "test",
    }


def test_RNA_seq_wrong_char(example_parameters):
    """Assert that RNA sequence contains only allowd characters"""
    sequences = {"five_end_flanking_sequence": "AGAACCGCACTCGGGAGGTTAAGCCGTCAATAACAATGACTCAAAGCCAAATGGTACTTGATAAGCGCAGCGTGGTGT",
                 "gene_of_interest": "bbbbbbbbbbbbbbbbbb",
                 "three_end_flanking_sequence": "CATCGCGGGAATGCTACCTAGACTTGGGGTAAGTCGATCCAGTGTCCC"}
    assert_wrong_char(example_parameters, sequences)


def test_five_end_seq_wrong_char(example_parameters):
    """Assert that five flanking sequence contains only allowd characters"""
    sequences = {"five_end_flanking_sequence": "wwwwwwwww!",
                 "gene_of_interest": "CAGATCCGCATATGTCCTACTTCTGTCGGTGAG",
                 "three_end_flanking_sequence": "CATCGCGGGAATGCTACCTAGACTTGGGGTAAGTCGATCCAGTGTCCC"}
    assert_wrong_char(example_parameters, sequences)


def test_three_end_seq_wrong_char(example_parameters):
    """Assert that three flanking sequence contains only allowd characters"""
    sequences = {"five_end_flanking_sequence": "AGAACCGCACTCGGGAGGTTAAGCCGTCAATAACAATGACTCAAAGCCAAATGGTACTTGATAAGCGCAGCGTGGTGT",
                 "gene_of_interest": "CAGATCCGCATATGTCCTACTTCTGTCGGTGAG",
                 "three_end_flanking_sequence": "klkklkK"}
    assert_wrong_char(example_parameters, sequences)
```
    
# tests/test_mfe_optimization.py
```python
import RNA
import pytest
from dnachisel import DnaOptimizationProblem
from common.objectives.MFE_optimization import MinimizeMFE


@pytest.fixture
def parameters():
    return {
        'input_sequence': 'CCGTCGCGGCAGGTTATTATACCTCATTCCTTGGAGACATACAACTATCAATGGGACTTGAGGTTAAGGTATTCCCGCATGAACGCGTGTACTGAAAATATGAAGGCGAGGGCGGAAGCTTTCATTAGCGAGCACCTACAACGTTAGAGTTGGTCGTGTCTTGCTATGCGTCCAGCACATCTGTAAGCCGGTATAAGGCCAGGGGCGGTACATATCGTACAGATCTAGTACATGTTGATAACTTTCATCTGTCGTAGGAAGGCGGAGCCGCCCCTGACGGACGTAGAAAGGGGAATGGGCACTGAGACCCAGTGAGCCCCTTTTGCGTTCTTGGCAAATACCTAGACCTTCTGGTCGTCCTATCGTAATATCTCCTGATACTCATGACAGCAGGATAGCAGCCTGCAACCTCCATGTACTTCGTTGGATTCTTTCCGAGTCTCGTGTGAGTAGATGCTTTGGGGAGTTACCTCTAACACATGGCTTGTTTATTCGTAATTCGACTCCCATGCTTGCTTTTAAACGTCTGTCAACATGAACATTCTGGTCGCACGACGATTAAGAAAGGGAACTTCGTGTTGATGTAGTAGGATATAGCAG',
        'five_end': 'CCGTCG',
        'entropy_window': 30}


def optimization(parameters, method):
    """Check that the optimized MFE is closer to zero than input MFE"""
    input_sequence = parameters['input_sequence']
    five_end = parameters['five_end']
    entropy_window = parameters['entropy_window']
    objective = MinimizeMFE(entropy_window=entropy_window, five_end=five_end, mfe_method=method)
    problem = DnaOptimizationProblem(input_sequence, objectives=[objective])
    problem.resolve_constraints()
    problem.optimize()

    optimized_sequence = problem.sequence
    _, input_mfe = RNA.fold(five_end + input_sequence[:entropy_window])
    _, output_mfe = RNA.fold(five_end + optimized_sequence[:entropy_window])

    assert input_mfe < output_mfe


def test_stem_loop_optimization(parameters):
    optimization(parameters, 'stem-loop')


def test_rna_fold_optimization(parameters):
    optimization(parameters, 'RNAfold')

```
    
# tests/test_optimization.py
```python
from Bio.Seq import Seq
from common.OptimizationProblems import initialize_optimization_problem
from common.OptimizationTask import optimization_task
from common.utils.Datatypes import OptimizationParameters
from common.Evaluation import Evaluation


def get_rare_codons(species, threshold):
    if threshold is None:
        return []
    # TODO: get tables for different species, if needed (e.g. from python_codon_tables)
    freqs = {"UAA": 0.3, "UAG": 0.24, "UGA": 0.47, "GCA": 0.23, "GCC": 0.4, "GCG": 0.11, "GCU": 0.27, "UGC": 0.54,
             "UGU": 0.46, "GAC": 0.54, "GAU": 0.46, "GAA": 0.42, "GAG": 0.58, "UUC": 0.54, "UUU": 0.46, "GGA": 0.25,
             "GGC": 0.34, "GGG": 0.25, "GGU": 0.16, "CAC": 0.58, "CAU": 0.42, "AUA": 0.17, "AUC": 0.47, "AUU": 0.36,
             "AAA": 0.43, "AAG": 0.57, "CUA": 0.07, "CUC": 0.2, "CUG": 0.4, "CUU": 0.13, "UUA": 0.08, "UUG": 0.13,
             "AUG": 1.0, "AAC": 0.63, "AAU": 0.47, "CCA": 0.28, "CCC": 0.32, "CCG": 0.11, "CCU": 0.29, "CAA": 0.27,
             "CAG": 0.73, "AGA": 0.21, "AGG": 0.21, "CGA": 0.11, "CGC": 0.18, "CGG": 0.2, "CGU": 0.08, "AGC": 0.24,
             "AGU": 0.15, "UCA": 0.15, "UCC": 0.22, "UCG": 0.05, "UCU": 0.19, "ACA": 0.28, "ACC": 0.36, "ACG": 0.11,
             "ACU": 0.25, "GUA": 0.12, "GUC": 0.24, "GUG": 0.46, "GUU": 0.18, "UGG": 1, "UAC": 0.56, "UAU": 0.44}
    rare_codons = [codon for codon, freq in freqs.items() if freq < threshold]
    return rare_codons


def test_optimization():
    """
    Test whole optimization process, check if optimized sequences conform to constraints set by parameters
    """
    seq1 = "AGCAGAGAAGGCGGAAGCAGTGGCGTCCGCAGCTGGGGCTTGGCCTGCGGGCGGCCAGCGAAGGTGGCGAAGGCTCCCACTGGATCCAGAGTTTGCCGTCCAAGCAGCCTCGTCTCGGCGCGCAGTGTCTGTGTCCGTCCTCTACCAGCGCCTTGGCTGAGCGGAGTCGTGCGGTTGGTGGGGGAGCCCTGCCCTCCTGGTTCGGCCTCCCCGCGCACTAGAACGATCATGAACTTCTGAAGGGACCCAGCTTTCTTTGTGTGCTCCAAGTGATTTGCACAAATAATAATATATATATTTATTGAAGGAGAGAATCAGAGCAAGTGATAATCAAGTTACTATGAGTCTGCTAAACTGTGAAAACAGCTGTGGATCCAGCCAGTCTGAAAGTGACTGCTGTGTGGCCATGGCCAGCTCCTGTAGCGCTGTAACAAAAGATGATAGTGTGGGTGGAACTGCCAGCACGGGGAACCTCTCCAGCTCATTTATGGAGGAGATCCAGGGATATGATGTAGAGTTTGACCCACCCCTGGAAAGCAAGTATGAATGCCCCATCTGCTTGATGGCATTACGAGAAGCAGTGCAAACGCCATGCGGCCATAGGTTCTGCAAAGCCTGCATCATAAAATCAATAAGGGATGCAGGTCACAAATGTCCAGTTGACAATGAAATACTGCTGGAAAATCAACTATTTCCAGACAATTTTGCAAAACGTGAGATTCTTTCTCTGATGGTGAAATGTCCAAATGAAGGTTGTTTGCACAAGATGGAACTGAGACATCTTGAGGATCATCAAGCACATTGTGAGTTTGCTCTTATGGATTGTCCCCAATGCCAGCGTCCCTTCCAAAAATTCCATATTAATATTCACATTCTGAAGGATTGTCCAAGGAGACAGGTTTCTTGTGACAACTGTGCTGCATCAATGGCATTTGAAGATAAAGAGATCCATGACCAGAACTGTCCTTTGGCAAATGTCATCTGTGAATACTGCAATACTATACTCATCAGAGAACAGATGCCTAATCATTATGATCTAGACTGCCCTACAGCCCCAATTCCATGCACATTCAGTACTTTTGGTTGCCATGAAAAGATGCAGAGGAATCACTTGGCACGCCACCTACAAGAGAACACCCAGTCACACATGAGAATGTTGGCCCAGGCTGTTCATAGTTTGAGCGTTATACCCGACTCTGGGTATATCTCAGAGGTCCGGAATTTCCAGGAAACTATTCACCAGTTAGAGGGTCGCCTTGTAAGACAAGACCATCAAATCCGGGAGCTGACTGCTAAAATGGAAACTCAGAGTATGTATGTAAGTGAGCT"
    seq2 = "GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGCCTTCACCCTCTGCTCTGGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTATCATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATCCTTCCTTGCAGGAAACCAGTCTCAGTGTCCAACTCTCTAACCTTGGAACTGTGAGAACTCTGAGGACAAAGCAGCGGATACAACCTCAAAAGACGTCTGTCTACATTGAATTGGGATCTGATTCTTCTGAAGATACCGTTAATAAGGCAACTTATTGCAGTGTGGGAGATCAAGAATTGTTACAAATCACCCCTCAAGGAACCAGGGATGAAATCAGTTTGGATTCTGCAAAAAAGGCTGCTTGTGAATTTTCTGAGACGGATGTAACAAATACTGAACATCATCAACCCAGTAATAATGATTTGAACACCACTGAGAAGCGTGCAGCTGAGAGGCATCCAGAAAAGTATCAGGGTAGTTCTGTTTCAAACTTGCATGTGGAGCCATGTGGCACAAATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACTAAAGACAGAATGAATGTAGAAAAGGCTGAATTCTGTAATAAAAGCAAACAGCCTGGCTTAGCAAGGAGCCAACATAACAGATGGGCTGGAAGTAAGGAAACATGTAATGATAGGCGGACTCCCAGCACAGAAAAAAAGGTAGATCTGAATGCTGATCCCCTGTGTGAGAGAAAAGAATGGAATAAGCAGAAACTGCCATGCTCAGAGAATCCTAGAGATACTGAAGATGTTCCTTGGATAACACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACATGATGGGGAGTCTGAATCAAATGCCAAAGTAGCTGATGTATTGGACGTTCTAAATGAGGTAGATGAATATTCTGGTTCTTCAGAGAAAATAGACTTACTGGCCAGTGATCCTCATGAGGCTTTAATATGTAAAAGTGAAAGAGTTCACTCCAAATCAGTAGAGAGTAATATTGAAGACAAAATATTTGGGAAAACCTATCGGAAGAAGGCAAGCCTCCCCAACTTAAGCCATGTAACTGAAAATCTAATTATAGGAGCATTTGTTACTGAGCCACAGATAATACAAGAGCGTCCCCTCACAAATAAATTAAAGCGTAAAAGGAGACCTACATCAGGCCTTCATCCTGAGGATTTTATCAAGAAAGCAGATTTGGCAGTTCAAAAGACTCCTGAAATGATAAATCAGGGAACTAACCAAACGGAGCAGAATGGTCAAGTGATGAATATTACTAATAGTGGTCATGAGAATAAAACAAAAGGTGATTCTATTCAGAATGAGAAAAATCCTAACCCAATAGAATCACTCGAAAAAGAATCTGCTTTCAAAACGAAAGCTGAACCTATAAGCA"
    seq3 = "GAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTAGAGTTCATCCTGGCAGCTCACCGACCTGCCTGCAAGATCCCTGATGACCTGGGCTTCCCAGAAGAGATGTCTGTGGCTTCCCTTGATCTGACTGGGGGCCTGCCAGAGGTTGCCACCCCGGAGTCTGAGGAGGCCTTCACCCTGCCTCTCCTCAATGACCCTGAGCCCAAGCCCTCAGTGGAACCTGTCAAGAGCATCAGCAGCATGGAGCTGAAGACCGAGCCCTTTGATGACTTCCTGTTCCCAGCATCATCCAGGCCCAGTGGCTCTGAGACAGCCCGCTCCGTGCCAGACATGGACCTATCTGGGTCCTTCTATGCAGCAGACTGGGAGCCTCTGCACAGTGGCTCCCTGGGGATGGGGCCCATGGCCACAGAGCTGGAGCCCCTGTGCACTCCGGTGGTCACCTGTACTCCCAGCTGCACTGCTTACACGTCTTCCTTCGTCTTCACCTACCCCGAGGCTGACTCCTTCCCCAGCTGTGCAGCTGCCCACCGCAAGGGCAGCAGCAGCAATGAGCCTTCCTCTGACTCGCTCAGCTCACCCACGCTGCTGGCCCTGTGAGGGGGCAGGGAAGGGGAGGCAGCCGGCACCCACAAGTGCCACTGCCCGAGCTGGTGCATTACAGAGAGGAGAAACACATCTTCCCTAGAGGGTTCCTGTAGACCTAGGGAGGACCTTATCTGTGCGTGAAACACACCAGGCTGTGGGCCTCAAGGACTTGAAAGCATCCATGTGTGGACTCAAGTCCTTACCTCTTCCGGAGATGTAGCAAAACGCATGGAGTGTGTATTGTTCCCAGTGACACTTCAGAGAGCTGGTAGTTAGTAGCATGTTGAGCCAGGCCTGGGTCTGTGTCTCTTTTCTCTTTCTCCTTAGTCTTCTCATAGCATTAACTAATCTATTGGGTTCATTATTGGAATTAACCTGGTGCTGGATATTTTCAAATTGTATCTAGTGCAGCTGATTTTAACAATAACTACTGTGTTCCTGGCAATAGTGTGTTCTGATTAGAAATGACCAATATTATACTAAGAAAAGATACGACTTTATTTTCTGGTAGATAGAAATAAATAGCTATATCCATGTACTGTAGTTTTTCTTCAACATCAATGTTCATTGTAATGTTACTGATCATGCATTGTTGAGGTGGTCTGAATGTTCTGACATTAACAGTTTTCCATGAAAACGTTTTATTGTGTTTTTAATTTATTTATTAAGATGGATTCTCAGATATTTATATTTTTATTTTATTTTTTTCTACCTTGAGGTCTTTTGACATGTGGAAAGTGAATTTGAATGAAAAATTTAAGCATTGTTTGCTTATTGTTCCAAGACATTGTCAATAAAAGCATTTAAGTTGAATGC"
    sequences = [seq1, seq2, seq3]

    for seq in sequences:
        parameters = OptimizationParameters(
            input_mRNA=seq,
            input_DNA=seq.replace('U', 'T'),
            five_end="ATG",
            three_end="CGG",
            avoid_codons=['ACT', 'CTA'],
            avoid_motifs=[],
            max_GC_content=0.7,
            min_GC_content=0.3,
            GC_window_size=30,
            usage_threshold=0.1,
            uridine_depletion=True,
            organism='h_sapiens',
            entropy_window=80,
            number_of_sequences=1,
            mfe_method='stem_loop',
            dinucleotides=True,
            codon_pair=False,
            CAI=False,
            location=(0, len(seq) - len(seq) % 3, 1),
            filename="optimization_results"
        )
        translation = Seq(parameters.five_end + seq + parameters.three_end).translate()
        optimization_problem = initialize_optimization_problem(parameters)
        sequence = optimization_task(optimization_problem)
        evaluator = Evaluation([sequence], parameters)
        seq_properties = evaluator.get_evaluation()
        for properties in seq_properties["optimized"]:
            opt_seq = properties["RNASeq"]
            u_depletion = sum(1 for i in range(0, len(opt_seq) - 2, 3) if opt_seq[i + 2] in ["U", "T"])
            avoid_motifs = sum(1 for motif in parameters.avoid_motifs if motif in opt_seq)
            assert parameters.min_GC_content <= properties["GC_ratio"] <= parameters.max_GC_content
            # assert properties['score'] > seq_properties['input']['score']
            # This assertion fails
            if parameters.uridine_depletion:
                assert u_depletion == 0
            rare_codons = get_rare_codons(parameters.organism, parameters.usage_threshold)
            breaches = sum(1 for i in range(0, len(opt_seq) - 2, 3) if opt_seq[i:i + 2] in rare_codons)
            assert Seq(opt_seq).translate() == translation
            assert avoid_motifs == 0
            assert breaches == 0

```
    
# tests/test_organism.py
```python
from common.utils.RequestParser import RequestParser
from common.utils.Exceptions import SpeciesError
import pytest
from tests.test_request_parsing import MockRequest


@pytest.fixture
def example_parameters():
    def _example_parameters(organism):
        return {
            "config": {
                "avoid_codons": [],
                "avoided_motifs": [],
                "codon_usage_frequency_threshold": 0.1,
                "max_GC_content": 0.7,
                "GC_window_size": 30,
                "min_GC_content": 0.3,
                "organism": organism,
                "entropy_window": 6,
                "number_of_sequences": 5
            },
            "dinucleotides": False,
            "match_codon_pair": False,
            "uridine_depletion": True,
            "faster_MFE_algorithm": False,
            "precise_MFE_algorithm": False,
            "CAI": False,
            "file_name": "test",
            "sequences": {
                "five_end_flanking_sequence": "AGAACCGCACTCGGGAGGTTAAGCCGTCAATAACAATGACTCAAAGCCAAATGGTACTTGATAAGCGCAGCGTGGTGT",
                "gene_of_interest": "ATACTTTTCTGCATCGACTAAGTCTGAACGGAGAACCCGGAAACATATGTTGCTTCAGAACGCCGAGATTACTGCGTTTAGACTCAAACGAGAAGCGCCCCGCTCGTGTTATCTTCCAACTACCACGTATGTGTTGGACTGTTAGTGCACTTCAGATTGGTCCAAAACACCCTAGATGCTGAACCTGTACCCAGTAAGAAGTACGTAGTCGTTACAAGCAAACCCTAGCTTAGTAAAGGAAGCGCGTCTTACGTGTATGGTGGTTATCTTCCACCGTTGGACTTACGCGATTTGTGGAAATTACTTGCTTGCGTTGCTACGGGATTGAGAAAGCTTGCATATTGGGGGGTCAGACCCGTCATACAGGCCATCATGTCGGAAGGGCAGGTTCCATTGCCAGGTAGCTCCACCGCCCTGACGGGTTTCCAGGCACAGGTATACTCTCATACCAACAAAGCGCTTTCTGGGACTCGACGACACCAGGCCACGCGATGCTTCTCTAGATTCTCTGCGTTCGCTGTTGCATGCAGCCTCGCCTCGTTGAAGATGAGTATACAAACCCTGAACTGGGTGGTCCGATTGCTTGTAGGGTACACAACATGAGGGCGCACAAGCCCGGCGGGCCAATCCTTGATCACAACATATTCTACCTTCACTCCGCCAGACTGCGACCGCGTCAGTGTAAGCCCAAGTTTCGTGTTAGAAATACCACATGAATCGCGATTAACGGAGCACACGGCGCACGACGTCGATCCCTACTTAAGGGGTGTTTGACTCTCACAGTATGTTTATGCTCCGTACAGGCGCCCTAGTCGCAAGAAAGCCACTATGCTTGCTAATTTTTTCTTTGGCACGAACGTGCTTTCTAACTGGTAGATCCTATACTAGGCCTACTCGCTA",
                "three_end_flanking_sequence": "CATCGCGGGAATGCTACCTAGACTTGGGGTAAGTCGATCCAGTGTCCC"
            }
        }

    return _example_parameters


def test_correct_organism_human(example_parameters):
    """Assert correct human organism name is given"""
    organism_var = ["h_sapiens", "sapiens", "human", "homo_sapiens", "h sapiens", "homo sapiens", "H sapiens", "Homo Sapiens", "H_sapiens"]
    for organism in organism_var:
        organism_parametrs = example_parameters(organism=organism)
        request = MockRequest(**organism_parametrs)
        parser = RequestParser(request)
        parameters = parser.parse()
        assert parameters.organism == "h_sapiens"


def test_correct_organism_mouse(example_parameters):
    """Assert correct mouse organism name is given"""
    organism_var = ["m_musculus", "mouse", "mus_musc", "mus_musculus", "mus", "mus musculus", "m musculus", "M_musculus", "Mus Musculus"]
    for organism in organism_var:
        organism_parametrs = example_parameters(organism=organism)
        request = MockRequest(**organism_parametrs)
        parser = RequestParser(request)
        parameters = parser.parse()
        assert parameters.organism == "m_musculus"


def test_species_error_input(example_parameters):
    """Assert SpeciesError is raised when invalid organism is given"""
    request = MockRequest(**example_parameters(organism="wrong_organism"))
    parser = RequestParser(request)
    try:
        parameters = parser.parse()
    except SpeciesError as e:
        right_error_raised = True
    else:
        right_error_raised = False
    assert right_error_raised

```
    
# tests/test_request_parsing.py
```python
from unittest import main, TestCase
from unittest.mock import MagicMock as Mock
from common.utils.RequestParser import RequestParser
from json import dumps
MockRequest = Mock()

# Configure MockRequest to have the attributes and methods you need
MockRequest.return_value.get_data.return_value = dumps({
    "config": {
        "avoid_codons": ["ATG"],
        "avoided_motifs": [],
        "codon_usage_frequency_threshold": 0.1,
        "max_GC_content": 0.7,
        "GC_window_size": 30,
        "min_GC_content": 0.3,
        "organism": 'h_sapiens',
        "entropy_window": 6,
        "number_of_sequences": 5
    },
    "dinucleotides": False,
    "uridine_depletion": False,
    "faster_MFE_algorithm": False,
    "file_name": "test",
    "sequences": {
        "five_end_flanking_sequence": "ACG",
        "gene_of_interest": "AGCAAAGAAGGCGGA",
        "three_end_flanking_sequence": "ACG"
    },
    "precise_MFE_algorithm": None,
    "optimization_criterion": 'cai'
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
```
    
# tests/test_score.py
```python
from common.Evaluation import Evaluation
from common.utils.Datatypes import OptimizationParameters
import pytest


@pytest.fixture
def example_parameters():
    def _example_parameters(mRNA_seq, entropy_window):
        dna = mRNA_seq.replace('U', 'T')
        parameters = OptimizationParameters(
            input_mRNA=mRNA_seq,
            input_DNA=dna,
            five_end="ATG",
            three_end="CGG",
            avoid_codons=[],
            avoid_motifs=[],
            max_GC_content=0.7,
            min_GC_content=0.3,
            GC_window_size=30,
            usage_threshold=0.1,
            uridine_depletion=True,
            organism='h_sapiens',
            entropy_window=entropy_window,
            number_of_sequences=3,
            mfe_method='stem_loop',
            dinucleotides=False,
            location=(0, len(mRNA_seq) - len(mRNA_seq) % 3, 1),
            filename="optimization_results"
        )
        return parameters
    return _example_parameters


def test_ranking(example_parameters):
    """Check if the sequences are ranked properly by the Evaluation class"""
    seq_in = "AUGCCUGCGGAUCCGAGAAACGCCGGAAAGGACAAAGGAAAUCUUCAGAUGGAGCUCCAGAACUUGCUCACUUCCGCAAGAAUAAGGAAGAUGGAGUACGAAGCAAUAGUGGAAGAGCUUGAGGAAGACGAGCUCAAAUACGACCUCUUUGAGUACCAAGACUACCUUGAAAAUUACAUUAUGCCCUACGUUGAAAGAGCUUACAAAAACGGCAGUAAAGAACUGAUAGGAAUGGCGGAGGAAGUAAAGAGUAUAUUUGAAGAAAUUAUUGAUAUGAUAAAGGCGAGAAUAGAAAAGAAGUGA"
    seq1 = "AUGCCGGCGGACCCAAGAAACGCCGGCAAGGACAAGGGAAACUUGCAGAUGGAGCUGCAGAACUUGUUGACCUCCGCCCGAAUCAGGAAGAUGGAAUACGAAGCAAUAGUAGAGGAGCUGGAGGAAGACGAGCUGAAAUACGACCUCUUCGAGUACCAAGACUACCUCGAAAACUACAUCAUGCCCUACGUGGAACGCGCCUACAAAAACGGCUCAAAAGAGCUGAUAGGGAUGGCGGAGGAGGUCAAGAGCAUCUUCGAAGAGAUCAUCGACAUGAUAAAGGCACGGAUAGAAAAGAAAUGA"
    seq2 = "AUGCCAGCAGACCCGCGGAACGCCGGCAAGGACAAAGGAAACCUGCAGAUGGAACUGCAGAACUUGCUCACCUCCGCACGAAUAAGGAAGAUGGAGUACGAGGCAAUAGUGGAGGAGCUGGAGGAAGACGAGUUGAAAUACGACCUCUUCGAGUACCAAGACUACCUGGAAAACUACAUCAUGCCCUACGUCGAAAGAGCGUACAAAAACGGCUCAAAAGAACUGAUAGGGAUGGCCGAGGAAGUAAAGAGCAUCUUCGAAGAGAUCAUCGACAUGAUCAAGGCCCGCAUAGAGAAGAAGUGA"
    seq3 = "AUGCCGGCAGACCCACGAAACGCCGGCAAGGACAAAGGGAACCUGCAGAUGGAGCUGCAGAACUUGCUGACCUCCGCACGCAUAAGGAAGAUGGAGUACGAAGCAAUCGUGGAAGAGCUGGAGGAAGACGAGUUGAAAUACGACCUCUUCGAGUACCAAGACUACCUCGAAAACUACAUCAUGCCCUACGUCGAAAGAGCCUACAAAAACGGCAGCAAAGAACUGAUCGGAAUGGCCGAGGAGGUAAAGUCAAUAUUCGAGGAGAUCAUCGACAUGAUAAAGGCGCGGAUAGAAAAGAAGUGA"
    optimized = [seq1, seq2, seq3]
    parameters = example_parameters(seq_in, 10)

    evaluator = Evaluation(optimized, parameters)

    seq0_properties = evaluator.get_seq_properties("seq0", seq_in)
    seq1_properties = evaluator.get_seq_properties("seq1", seq1)
    seq2_properties = evaluator.get_seq_properties("seq2", seq2)
    seq3_properties = evaluator.get_seq_properties("seq3", seq3)

    # Sequence 1 is the best one
    assert seq0_properties.score < seq1_properties.score
    assert seq0_properties.score < seq2_properties.score
    assert seq0_properties.score < seq3_properties.score
    assert seq1_properties.score > seq2_properties.score
    assert seq2_properties.score > seq3_properties.score


def test_ranking_w_longer_window(example_parameters):
    """Check if the sequences are ranked properly by the Evaluation class"""
    seq_in = "AUGCCUGCGGAUCCGAGAAACGCCGGAAAGGACAAAGGAAAUCUUCAGAUGGAGCUCCAGAACUUGCUCACUUCCGCAAGAAUAAGGAAGAUGGAGUACGAAGCAAUAGUGGAAGAGCUUGAGGAAGACGAGCUCAAAUACGACCUCUUUGAGUACCAAGACUACCUUGAAAAUUACAUUAUGCCCUACGUUGAAAGAGCUUACAAAAACGGCAGUAAAGAACUGAUAGGAAUGGCGGAGGAAGUAAAGAGUAUAUUUGAAGAAAUUAUUGAUAUGAUAAAGGCGAGAAUAGAAAAGAAGUGA"
    seq1 = "AUGCCAGCAGACCCAAGAAAUGCAGGAAAAGACAAAGGAAAUCUGCAGAUGGAGCUCCAGAACUUGCUUACCUCCGCUCGGAUAAGGAAGAUGGAGUACGAAGCGAUCGUGGAGGAGCUUGAGGAAGACGAGCUGAAGUAUGAUCUCUUUGAGUACCAAGAUUACCUGGAGAACUAUAUUAUGCCCUACGUGGAACGAGCUUAUAAAAACGGCAGCAAGGAACUGAUAGGGAUGGCCGAGGAGGUCAAGUCUAUCUUCGAAGAAAUCAUUGAUAUGAUUAAGGCCCGCAUCGAAAAGAAAUGA"
    seq2 = "AUGCCAGCAGACCCAAGAAAUGCAGGAAAAGACAAAGGAAAUUUGCAGAUGGAGCUGCAGAACCUGCUCACCUCCGCCCGGAUAAGGAAGAUGGAGUACGAAGCUAUCGUGGAAGAGCUUGAGGAAGAUGAACUGAAAUACGACCUCUUCGAAUAUCAAGAUUACCUCGAGAACUAUAUUAUGCCUUACGUCGAGCGCGCGUAUAAGAACGGCUCUAAAGAGCUGAUCGGGAUGGCUGAAGAAGUUAAGAGCAUCUUUGAGGAGAUUAUUGAUAUGAUCAAGGCCCGAAUAGAGAAGAAGUGA"
    seq3 = "AUGCCAGCAGACCCAAGAAAUGCAGGAAAAGACAAAGGAAAUCUCCAGAUGGAGCUGCAGAACUUGCUUACCUCCGCCCGAAUUAGGAAGAUGGAGUACGAAGCCAUCGUAGAAGAGCUUGAGGAGGACGAGCUCAAGUAUGAUCUGUUUGAAUACCAAGAUUAUCUGGAGAACUAUAUUAUGCCCUACGUUGAACGGGCUUACAAAAACGGCUCUAAGGAACUGAUCGGGAUGGCUGAGGAGGUCAAAAGCAUAUUCGAAGAGAUUAUCGAUAUGAUCAAGGCGCGCAUAGAAAAGAAGUGA"
    optimized = [seq1, seq2, seq3]
    parameters = example_parameters(seq_in, 40)

    evaluator = Evaluation(optimized, parameters)

    seq0_properties = evaluator.get_seq_properties("seq0", seq_in)
    seq1_properties = evaluator.get_seq_properties("seq1", seq1)
    seq2_properties = evaluator.get_seq_properties("seq2", seq2)
    seq3_properties = evaluator.get_seq_properties("seq3", seq3)

    # Sequence 1 is the best one
    assert seq0_properties.score < seq1_properties.score
    assert seq0_properties.score < seq2_properties.score
    assert seq0_properties.score < seq3_properties.score
    assert seq1_properties.score > seq2_properties.score
    assert seq2_properties.score > seq3_properties.score


def test_gc_score(example_parameters):
    """Check if the sequences with optimal GC content are ranked higher than sequences with extreme GC content"""
    seq_gc = "GCCGGCCCCCGGTGCGACGAGGTGAGCCCCTGGAGCCGG"
    seq_at = "TGATTCATCAAAATGAACACTTACATCACTAAATGAAAA"
    seq_opt = "CGGTGAGACAAACTGAACTGGTTCGGCTACGAGAAACTG"
    seqs = [seq_gc, seq_opt]
    parameters = example_parameters(seq_at, 5)

    evaluator = Evaluation(seqs, parameters)

    seq_opt_properties = evaluator.get_seq_properties("seq_opt", seq_opt)
    seq_gc_properties = evaluator.get_seq_properties("seq_gc", seq_gc)
    seq_at_properties = evaluator.get_seq_properties("seq_at", seq_at)

    assert seq_opt_properties.score > seq_at_properties.score
    assert seq_opt_properties.score > seq_gc_properties.score



```
    
# tests/test_uridine_depletion.py
```python
from dnachisel import DnaOptimizationProblem
from dnachisel.Location import Location
from common.constraints.UridineDepletion import UridineDepletion


def test_uridine_depletion_with_occurrences():
    """Check if breach locations are correctly identified"""
    uridine_depletion = UridineDepletion()
    sequence = "AATCTAGTCCCTACGTGA"
    problem = DnaOptimizationProblem(sequence, [uridine_depletion])
    evaluation = uridine_depletion.evaluate(problem)
    assert len(evaluation.locations) == 2
    assert str(evaluation.locations[0]) == str(Location(0, 3, 1))
    assert str(evaluation.locations[1]) == str(Location(9, 12, 1))


def test_uridine_depletion_no_occurence():
    """Check if breach locations are correctly identified (none in this case)"""
    uridine_depletion = UridineDepletion()
    sequence = "AAGCTAGTCCCAACGTGA"
    problem = DnaOptimizationProblem(sequence, [uridine_depletion])
    evaluation = uridine_depletion.evaluate(problem)
    assert len(evaluation.locations) == 0
    assert evaluation.score == 0


def test_uridine_depletion_optimization_results():
    """Check if the uridine depletion is correctly applied during optimization (no U at third position in any codon)"""
    # optimization
    uridine_depletion = UridineDepletion()
    sequence = "AATCTAGTCCCTACGTGA"
    problem = DnaOptimizationProblem(sequence, [uridine_depletion])
    problem.resolve_constraints()
    problem.optimize()
    optimized_sequence = problem.sequence

    # evaluation
    for_evaluation = DnaOptimizationProblem(optimized_sequence, [uridine_depletion])
    evaluation = uridine_depletion.evaluate(for_evaluation)
    for i in range(0, len(optimized_sequence) - 2, 3):
        assert optimized_sequence[i + 2] != "T"
    assert len(evaluation.locations) == 0
    assert evaluation.score == 0

```
    
