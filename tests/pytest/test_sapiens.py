import numpy as np
import sapiens


def test_sapiens_predict():
    seq = 'QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNRVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS'
    pred = sapiens.predict_scores(
        seq=seq,
        chain_type='H',
        model_version='latest'
    )
    assert pred.shape == (len(seq), 20), 'Expected matrix (length of sequence * 20 amino acids)'
    assert (pred.idxmax(axis=1).values == np.array(list(seq))).sum() > 100, 'Prediction should be similar to input sequence'
