"""Removes probes with long stretches of A or T bases.

This acts as a filter on the probes by returning a list of
probes, among the input probes, that do _not_ have a long
stretch of A or T, tolerating some given number of
mismatches. It preserves the order of the input.
"""

from collections import OrderedDict

from catch.filter.base_filter import BaseFilter
from catch import probe
from catch.utils import longest_common_substring as lcf

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class PolyGFilter(BaseFilter):
    """Filter that removes probes with poly(G).
    """

    def __init__(self, length, mismatches, min_exact_length_to_consider=6):
        """
        Args:
            length/mismatches: remove probes that contain at least
                LENGTH 'G' bases in a row, tolerating up to MISMATCHES
                mismatches
            min_exact_length_to_consider: only look for a stretch of 'G'
                (according to length/mismatches) in a probe if that
                probe contains an exact stretch of 'G' that is
                at least this length long. This is only meant to improve
                runtime, because the call to
                probe.Probe.longest_common_substring_length() is slow;
                this can reduce the number of times that needs to be called,
                but also result in false negatives. To always look
                for a stretch according only to length/mismatches, set
                the value of this argument to 0.
        """
        self.length = length
        self.mismatches = mismatches
        self.min_exact_length_to_consider = min_exact_length_to_consider

    def _filter(self, input):
        """Return a subset of the input probes.
        """
        if len(input) == 0:
            return input

        exact_g_stretch = 'G'*self.min_exact_length_to_consider

        probe_len = len(input[0])
        for p in input:
            probe_len = max(probe_len, len(p))
        g_stretch = probe.Probe.from_str('G'*probe_len)

        out = []
        for p in input:
            keep = True
            if exact_g_stretch in p.seq_str:
                lcf_len = p.longest_common_substring_length(
                    g_stretch, self.mismatches)
                if lcf_len >= self.length:
                    # The stretch exceeds the limit
                    keep = False
            if keep:
                out += [p]
        return out
