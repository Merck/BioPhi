import math
import torch
from fairseq import metrics, modules, utils
from fairseq.criterions import register_criterion, FairseqCriterion
from fairseq.criterions.label_smoothed_cross_entropy import label_smoothed_nll_loss

@register_criterion('smooth_masked_lm')
class SmoothMaskedLmLoss(FairseqCriterion):
    """
    Label-smoothed implementation for the loss used in masked language model (MLM) training.
    """

    def __init__(self, task, label_smoothing):
        super().__init__(task)
        self.eps = label_smoothing

    @staticmethod
    def add_args(parser):
        """Add criterion-specific arguments to the parser."""
        # fmt: off
        parser.add_argument('--label-smoothing', default=0., type=float, metavar='D',
                            help='epsilon for label smoothing, 0 means no label smoothing')
        # fmt: on

    def forward(self, model, sample, reduce=True):
        """Compute the loss for the given sample.

        Returns a tuple with three elements:
        1) the loss
        2) the sample size, which is used as the denominator for the gradient
        3) logging outputs to display while training
        """
        # compute MLM loss
        masked_tokens = sample['target'].ne(self.padding_idx)

        # Rare: when all tokens are masked, project all tokens.
        # We use torch.where to avoid device-to-host transfers,
        # except on CPU where torch.where is not well supported
        # (see github.com/pytorch/pytorch/issues/26247).
        if masked_tokens.device == torch.device('cpu'):
            if not masked_tokens.any():
                masked_tokens.fill_(True)
        else:
            masked_tokens = torch.where(
                masked_tokens.any(),
                masked_tokens,
                masked_tokens.new([True]),
            )

        net_output = model(**sample['net_input'], masked_tokens=masked_tokens)
        logits, _ = net_output
        targets = model.get_targets(sample, [logits])
        targets = targets[masked_tokens]

        lprobs = model.get_normalized_probs(net_output, log_probs=True)

        loss, nll_loss = label_smoothed_nll_loss(
            lprobs.view(-1, lprobs.size(-1)),
            targets.view(-1, 1),
            self.eps,
            reduce=reduce,
            ignore_index=self.padding_idx
        )

        sample_size = masked_tokens.int().sum()
        logging_output = {
            'loss': loss.data,
            'nll_loss': nll_loss.data,
            'ntokens': sample['ntokens'],
            'nsentences': sample['nsentences'],
            'sample_size': sample_size,
        }
        return loss, sample_size, logging_output

    @staticmethod
    def reduce_metrics(logging_outputs) -> None:
        """Aggregate logging outputs from data parallel training."""
        loss_sum = sum(log.get('loss', 0) for log in logging_outputs)
        nll_loss_sum = sum(log.get('nll_loss', 0) for log in logging_outputs)
        sample_size = sum(log.get('sample_size', 0) for log in logging_outputs)

        metrics.log_scalar('loss', loss_sum / sample_size / math.log(2), sample_size, round=3)
        metrics.log_scalar('nll_loss', nll_loss_sum / sample_size / math.log(2), sample_size, round=3)
        metrics.log_derived('ppl', lambda meters: utils.get_perplexity(meters['nll_loss'].avg))

    @staticmethod
    def logging_outputs_can_be_summed() -> bool:
        """
        Whether the logging outputs returned by `forward` can be summed
        across workers prior to calling `reduce_metrics`. Setting this
        to True will improves distributed training speed.
        """
        return True
