Alan ran the first 50 simulations in the parameter file, and a number of them didn't finish within 2 weeks. 
I hypothesize that this issue arose due to failure of the pyslim script to coalesce, which has happened in the past when the migration matrix was not
properly specified.

In the:

20210818_150_V2_SummaryReport.html put together by Alan, all of the demographies that failed were Estuary.

On Discovery I looked at the outputs in:

`MVP-NonClinalAF/sim-output_150_V2`

'ls -lah`
revealed an interesting pattern - the dates the files are written seem to be in succession. Why is that? If the partition wasn’t busy, shouldn’t they have all started at the same time? (maybe this is why some of them timed out - the 2 weeks is from the time the array script was submitted)

Thoughts:
- 1. I need to check if the VCF files that I'm producing are too large to work with, even after filtering. Maybe reduce the population size or mutation rate in pyslim.
  - try opening VCF files in VCFR - there are some edits I wanted to make anyway
- 2. check that pyslim can work on the failed sims to produce a VCF file.
  - 1231143, 1231141, 1231096, 1231098 as starting points

