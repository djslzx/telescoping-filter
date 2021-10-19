# Telescoping Adaptive Filter (TAF)
David J. Lee, Samuel McCauley, Shikha Singh, Max Stein

This work appeared at ESA 2021.  If you use this software, please cite us:
```
@InProceedings{lee_et_al:LIPIcs.ESA.2021.60,
  author =	{Lee, David J. and McCauley, Samuel and Singh, Shikha and Stein, Max},
  title =	{{Telescoping Filter: A Practical Adaptive Filter}},
  booktitle =	{29th Annual European Symposium on Algorithms (ESA 2021)},
  pages =	{60:1--60:18},
  series =	{Leibniz International Proceedings in Informatics (LIPIcs)},
  ISBN =	{978-3-95977-204-4},
  ISSN =	{1868-8969},
  year =	{2021},
  volume =	{204},
  editor =	{Mutzel, Petra and Pagh, Rasmus and Herman, Grzegorz},
  publisher =	{Schloss Dagstuhl -- Leibniz-Zentrum f{\"u}r Informatik},
  address =	{Dagstuhl, Germany},
  URL =		{https://drops.dagstuhl.de/opus/volltexte/2021/14641},
  URN =		{urn:nbn:de:0030-drops-146410},
  doi =		{10.4230/LIPIcs.ESA.2021.60},
  annote =	{Keywords: Filters, approximate-membership query data structures (AMQs), Bloom filters, quotient filters, cuckoo filters, adaptivity, succinct data structures}
}
```

## Overview
The Telescoping Adaptive Filter (TAF) is a practical, provably adaptive filter.  The TAF supports insertions and queries at high throughputs.  In addition, it sustains its false positive guarantees over sequences of insertions by fixing false positives as they occur.

For theoretical analysis, performance benchmarks, and architectural details, please see our paper:

["Telescoping Filter: A Practical Adaptive Filter"](https://arxiv.org/abs/2107.02866), European Symposium on Algorithms (ESA) 2021 by David J. Lee, Samuel McCauley, Shikha Singh, and Max Stein.

## API
The TAF supports the following operations:
- `taf_lookup(filter, elt)`: Returns whether `elt` is in the `filter`. Note that lookups may return false positives (a characteristic of all filters).
- `taf_insert(filter, elt)`: Insert `elt` into the `filter`. 
- `taf_clear(filter)`: Remove all elements from the `filter`.

### Basic usage
Here is a simple C example that instantiates a TAF, inserts an element, queries the element, and then deallocates the TAF:

```C
TAF* filter = new_taf(64);      // Create a filter with 64 empty slots
taf_insert(filter, 1);          // Insert 1
assert(taf_lookup(filter, 1));  // Query 1
taf_destroy(filter);            // Deallocate the filter
```

### More usage examples
To see more extensive usage examples, see the TAF's testing code in `taf.c`, following the macro `#ifndef TEST_TAF`.

### Other filters
We also provide three additional filter implementations.

#### Extension adaptive filter (exAF)
`exaf.*` contains the exAF, a practical implementation of Bender et al.'s [Broom filter](https://arxiv.org/abs/1711.01616) that leverages the TAF's core architecture. Its API mirrors the TAF's:
- `exaf_lookup(filter, elt)`
- `exaf_insert(filter, elt)`
- `exaf_clear(filter)`

#### Uncompressed TAF (uTAF)
`utaf.*` contains the uTAF, a variant of the TAF that does not use compression when fixing false positives. Its API also mirrors the TAF's:
- `utaf_lookup(filter, elt)`
- `utaf_insert(filter, elt)`
- `utaf_clear(filter)`

#### Rank-and-select quotient filter (RSQF)
`rsqf.*` contains a from-scratch implementation of Pandey et al.'s RSQF, the quotient filter architecture that undergirds the [Counting Quotient Filter (CQF)](https://github.com/splatlab/cqf).  The TAF, uTAF, and exAF are built using this RSQF implementation.
- `rsqf_lookup(filter, elt)`
- `rsqf_insert(filter, elt)`
- `rsqf_clear(filter)`

## Build and test
To build the TAF from `src/`:
```
make taf
```
To run the tests in `taf.c`, after building:
```
./taf.o
```

Similar `make` commands are available for `utaf`, `exaf`, `rsqf`, and `arcd`.

## Authors
- David J. Lee <djl328@cornell.edu>
- Samuel McCauley
- Shikha Singh
- Max Stein
