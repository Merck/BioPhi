import os

from abnumber import Chain

from biophi.humanization.methods.sapiens.roberta import RoBERTaSeq2Seq
# register loss function
# noinspection PyUnresolvedReferences
from biophi.humanization.methods.sapiens.smooth_masked_lm import SmoothMaskedLmLoss

MODEL_VERSIONS = ['v1']

CACHED_MODELS = {}


def load_model(model_dir, checkpoint_name, cache=True):
    model_key = (model_dir, checkpoint_name)
    if not cache or model_key not in CACHED_MODELS:
        CACHED_MODELS[model_key] = RoBERTaSeq2Seq.load(model_dir, checkpoint_name)
    return CACHED_MODELS[model_key]


def get_model_path(model_version):
    src_root = os.path.dirname(__file__)
    return os.path.join(src_root, 'models', model_version)


def sapiens_predict_chain(chain: Chain, model_version='latest', return_all_hiddens=False):
    return sapiens_predict_seq(
        seq=chain.seq,
        chain_type=chain.chain_type,
        model_version=model_version,
        return_all_hiddens=return_all_hiddens
    )


def sapiens_predict_seq(seq, chain_type, model_version='latest', return_all_hiddens=False):
    assert '/' not in model_version
    if model_version == 'latest':
        model_version = MODEL_VERSIONS[-1]
    model_dir = get_model_path(model_version)
    if chain_type == 'H':
        checkpoint_name = 'checkpoint_vh.pt'
    elif chain_type in ['K', 'L']:
        checkpoint_name = 'checkpoint_vl.pt'
    else:
        raise ValueError(f'Unknown chain type {chain_type}')
    model = load_model(model_dir, checkpoint_name)
    return model.predict_proba(seq, return_all_hiddens=return_all_hiddens)