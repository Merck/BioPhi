import numpy as np
import sapiens
from abnumber import Chain
from biophi.humanization.methods.humanization import humanize_chain, SapiensHumanizationParams


def test_sapiens_predict():
    seq = 'QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNRVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS'
    pred = sapiens.predict_scores(
        seq=seq,
        chain_type='H',
        model_version='latest'
    )
    assert pred.shape == (len(seq), 20), 'Expected matrix (length of sequence * 20 amino acids)'
    assert (pred.idxmax(axis=1).values == np.array(list(seq))).sum() > 100, 'Prediction should be similar to input sequence'


def test_sapiens_humanize():
    seq = 'QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNRVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS'
    chain = Chain(seq, scheme='kabat')
    humanization_params = SapiensHumanizationParams(
        model_version='latest',
        humanize_cdrs=False,
        scheme='kabat',
        cdr_definition='kabat',
        iterations=1
    )
    humanization = humanize_chain(chain, params=humanization_params)
    assert humanization.humanized_chain.seq == 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNRVTLTTDTSTTTAYMELRSLRSDDTAVYYCARRDYRFDMGFDYWGQGTLVTVSS'