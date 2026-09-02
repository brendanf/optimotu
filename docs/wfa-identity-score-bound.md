# Identity-feasible WFA band and score bounds

WFA clustering and `seq_distmx()` measure pairwise **CIGAR identity**, not the
affine alignment score:

\[
d = (L - M) / L
\]

where \(M\) is the number of match columns and \(L\) is the alignment length
(matches + mismatches + insertions + deletions). A pair is kept when
\(d \le \tau\). Constrained WFA still needs a **diagonal band** and a
**`setMaxAlignmentSteps` score cap**. Those limits must be *sound*: no
alignment with \(d \le \tau\) may require a wider band or a higher score.
This note derives those limits. Implementation:
`src/wfa_identity_bound.h`.

## Global accounting

Write \(s_1\) for the shorter sequence (\(l_1 \le l_2\)) and \(\sigma = 1 -
\tau\). Any global CIGAR satisfies

\[
M + X + D = l_1, \qquad M + X + I = l_2,
\]

so \(I = D + \Delta\) with \(\Delta = l_2 - l_1\), alignment length
\(L = l_2 + D\), and total non-match columns \(E = X + I + D\).

\(d \le \tau\) is equivalent to \(M \ge \sigma L\). Substituting the length
identities yields the length filter \(l_1 / l_2 \ge \sigma\) (already
applied before banding) and

\[
D \le \frac{l_1 - \sigma l_2}{1+\sigma}, \qquad
I \le \frac{l_2 - \sigma l_1}{1+\sigma}, \qquad
E \le \frac{\tau(l_1 + l_2)}{1+\sigma}.
\]

The last quantity is the historical `maxd1` used as an **edit** cap. It is
the maximum number of non-match columns, attained when \(X = 0\) and the
gap counts are maximal. It is *not* a valid affine score.

## Band (\(k\))

WFA’s diagonal \(k\) is a prefix of \(I - D\). A feasible CIGAR has at most
\(I_{\max}\) insertions and \(D_{\max}\) deletions in total, so every prefix
obeys

\[
-D_{\max} \le k \le I_{\max}.
\]

Those are the existing `min_k` / `max_k` formulas. They depend only on
identity and sequence lengths, not on gap penalties. Extend-span searches
use a wider negative band (`min_k = -\lceil l_1 \sigma \rceil`) because
leading gaps in the shorter sequence are scored but trailing end gaps are
stripped from \(d\).

## Score cap

`setMaxAlignmentSteps(S)` aborts when the WFA **score** is \(\ge S\). A
false abort happens only if some CIGAR with \(d \le \tau\) has score
\(\ge S\). The sound cap is therefore one more than

\[
S_{\max} = \max \{\mathrm{score}(A) : d(A) \le \tau\}
\]

over CIGAR structure (not sequence content). This is the affine analogue of
`maxd1`.

With WFA `match \le 0` (zero or a match bonus), extra matches cannot
increase score, so \(S_{\max}\) is attained on the identity boundary. Gap
characters of total length \(G\) are most expensive when they are split
into length-1 runs, or (for dual affine) when each run uses the more
costly of the two pieces. That worst-case gap cost is

\[
c_{\max}(G) = \max\bigl(G \cdot \max_p(o_p + e_p),\;
\max_p(o_p + e_p G)\bigr).
\]

The identity polytope in \((X, D)\) is linear. \(S\) is linear in those
counts after substituting \(c_{\max}\), so the maximum is at a vertex:

1. **Mismatch-heavy:** \(D = 0\), \(I = \Delta\),
   \(X_{\max} = \tau l_2 - \Delta\)
   \[
   S_1 = x\,X_{\max} + c_{\max}(\Delta)
   \]
2. **Gap-heavy:** \(X = 0\), \(G = E_{\max} = \texttt{maxd1}\)
   \[
   S_2 = c_{\max}(E_{\max})
   \]

Then \(S_{\max} = \max(S_1, S_2)\). Callers pass \(S_{\max} + 1\) because
WFA aborts on `score >= max_alignment_steps`.

For edit penalties \((x, o, e) = (1, 0, 1)\), both vertices collapse to
\(E_{\max}\), recovering `maxd1`.

If `match > 0` (unusual for WFA), matches are a per-column cost and
\(S_{\max}\) is increased by `match * l1`, which is a sound over-estimate.

## What the bound does not claim

WFA returns a **minimum-score** CIGAR, which is then converted to identity.
Affine-optimal \(\neq\) identity-optimal, so a completed alignment can still
have \(d > \tau\) even when some other CIGAR would pass. The cap only
guarantees that WFA is not aborted on a CIGAR that would have satisfied the
threshold. The length-1-run construction is pessimistic on purpose: it does
not use sequence content.

## Call sites

Constrained WFA in `Wfa2ClusterWorker`, `Wfa2DistWorker`,
`Wfa2SearchWorker`, `HybridClusterWorker`, `HybridSearchWorker`,
`kmer.cpp`, and `prealign_dist_matrix.cpp` uses `wfa_identity_bound()`.
Edlib’s `k` remains an edit-distance limit (`maxd1`). The edit prealigner
in `seq_distmx_prealign()` also still uses `maxd1`; only the final WFA
aligner uses the affine cap.
