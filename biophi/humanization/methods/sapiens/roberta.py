from fairseq.models import register_model_architecture
from fairseq.models.roberta import RobertaModel, RobertaHubInterface
from fairseq.models.roberta.model import base_architecture
from fairseq.tasks.sentence_prediction import SentencePredictionTask
from fairseq.tasks.masked_lm import MaskedLMTask
import fairseq
import fairseq.utils
import torch
import pandas as pd
from dataclasses import dataclass
import os


def seq2tokens(seq, dictionary, add_bos):
    assert ' ' not in seq, f'Expected regular sequence without space separators, found space in "{seq}"'
    sentence = ('<s> ' if add_bos else '') + seq2sentence(seq).replace('X', '<mask>') + ' </s>'
    return dictionary.encode_line(sentence, append_eos=False, add_if_not_exist=False).long()


def seq2sentence(seq):
    """Convert sequence to string of space-separated letters"""
    return ' '.join(list(seq))


@dataclass
class RoBERTa:
    interface: RobertaHubInterface

    @classmethod
    def load(cls, model_dir, checkpoint_name='checkpoint.pt') -> 'RoBERTa':
        full_model_dir = os.path.abspath(model_dir)
        if os.path.exists(full_model_dir):
            interface = RobertaModel.from_pretrained(
                full_model_dir,
                checkpoint_name,
                full_model_dir,
                bpe=None
            )
            interface.eval()
            return cls(interface=interface)
        else:
            raise FileNotFoundError(f'Model not found: {full_model_dir}')

    def predict_proba(self, seq):
        raise NotImplementedError()


@dataclass
class RoBERTaSeq2Seq(RoBERTa):

    def predict_proba(self, seq, remove_special=True, return_all_hiddens=False):
        """
        Get model output probabilities for all positions in a single sequence (slow, consider batching for large inputs)
        :param seq: Input sequence (str or Bio.Seq)
        :param remove_special: remove special dictionary tokens from predicted probability matrix
        :param return_all_hiddens: return attention from all layers (used with return_attention)
        :return: 2D dataframe of output token probabilities (columns) for each input position (rows)
        """
        add_bos = self._is_adding_bos()
        tokens = self._source_seq_tokens(seq, add_bos=add_bos)
        with torch.no_grad():
            logits, extra = self.interface.model(
                tokens.unsqueeze(0),
                features_only=False,
                return_all_hiddens=return_all_hiddens
            )
            pred = logits.softmax(dim=-1)[0]
            # remove EOS token
            pred = pred[:-1, :]
            if add_bos:
                # remove BOS token
                pred = pred[1:, :]
        pred = pd.DataFrame(pred.numpy(), columns=self.interface.task.target_dictionary.symbols)
        if remove_special:
            pred.drop(['<s>', '<pad>', '</s>', '<unk>', '<mask>'], axis=1, inplace=True)
        return pred

    def _is_adding_bos(self):
        if isinstance(self.interface.task, SentencePredictionTask):
            return True
        elif isinstance(self.interface.task, MaskedLMTask):
            return True
        return ValueError(f'Unknown task: {type(self.interface.task)}')

    def _source_seq_tokens(self, seq, add_bos):
        return seq2tokens(seq, dictionary=self.interface.task.source_dictionary, add_bos=add_bos)

    def _target_seq_tokens(self, seq, add_bos):
        return seq2tokens(seq, dictionary=self.interface.task.target_dictionary, add_bos=add_bos)


@register_model_architecture('roberta', 'roberta_small')
def roberta_small_architecture(args):
    args.encoder_layers = getattr(args, 'encoder_layers', 4)
    args.encoder_embed_dim = getattr(args, 'encoder_embed_dim', 128)
    args.encoder_ffn_embed_dim = getattr(args, 'encoder_ffn_embed_dim', 256)
    args.encoder_attention_heads = getattr(args, 'encoder_attention_heads', 8)
    base_architecture(args)
