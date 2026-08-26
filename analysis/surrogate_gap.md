# Surrogate gap

50 designed backbones, delta 1.0 nats.

The parser optimises the summed unconditional single-pass marginal. The model's
own likelihood is autoregressive and order dependent. This measures the distance
between them.

Conditional scores come from ProteinMPNN `--score_only`, mean negative
log-likelihood per residue under a random decoding order, averaged over ten
independent orders and multiplied by chain length to give totals in nats.

## Report the scatter, not the correlation

The correlation between the two scores is 0.735 across the whole shell and 0.060
inside the parser's top 500. Those are not two properties. They are one
relationship seen through two window widths, and a rank correlation inside a
narrow window is not interpretable on its own.

Applying Thorndike's case II correction forward, from the wide-window
relationship and the two spreads, predicts a within-band correlation of 0.033.
The observed 0.060 is slightly above that on 66 percent of backbones. Range
restriction accounts for the collapse and marginally over-explains it: the
surrogate retains a little more local signal than restriction alone would give,
not less.

The quantity that does not depend on window width is the residual scatter of the
conditional score about its regression on the surrogate.

| | |
|---|---|
| residual scatter | 2.42 nats (10th 1.98, 90th 2.70) |
| width of the parser's top-500 band | 0.30 nats |
| ratio | 7.9x |

**The conditional and unconditional scores differ by 2.4 nats root mean square,
roughly eight times the width of the parser's entire top-500 band. Ordering
within that band therefore carries no usable information about the conditional
model.** That statement is immune to the window width and is what should be
reported.

It also makes the earlier observation unsurprising rather than notable: there is
no k at which a larger candidate set stops helping, because widening k widens the
band, so the within-band correlation must rise. Measured, it climbs from -0.115
at k = 10 to 0.060 at k = 500, tracking the restriction prediction in shape.

## A correction to our own earlier reading

We previously argued that decoding-order noise could not explain the collapse,
and used that as the defence. That argument is about variance in the conditional
score and is a different question from the width of the window being correlated
inside. It is true that order noise is small, 0.023 per-residue nats against a
wide-window conditional range of 0.63, but it was never the relevant confound and
should not have been presented as the answer.

## What follows

Conditional rescoring is load-bearing rather than a refinement, and k should be
set by what can be afforded rather than by a recall target. The exactness
language survives with its qualifier: the parser is exact with respect to the
unconditional marginal, and the distance from that surrogate to the conditional
model is now reported as a scatter in nats rather than as a correlation.
