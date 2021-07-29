import numpy as np
from biophi.humanization.methods.sapiens.predict import sapiens_predict_seq


def test_sapiens_predict():
    seq = 'QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNRVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS'
    pred = sapiens_predict_seq(
        seq=seq,
        chain_type='H',
        model_version='latest'
    )
    assert pred.shape == (len(seq), 20), 'Expected matrix (length of sequence * 20 amino acids)'
    assert (pred.idxmax(axis=1).values == np.array(list(seq))).sum() > 100, 'Prediction should be similar to input sequence'
