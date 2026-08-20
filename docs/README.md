# cdr_FS documentation

The [README](../README.md) covers the method, the stages and an example run. The reference
material lives here, in six pages.

| Page | What is on it |
|---|---|
| [Quickstart](quickstart.md) | An annotated first run: what `check` tells you before you start, what the stage reports mean, where their numbers come from, and three outputs that look wrong and are not |
| [Troubleshooting](troubleshooting.md) | Exit codes and which stages refuse, what each label in a run's summary means, and where to look when a run retains nothing or a figure does not draw |
| [Configuration](configuration.md) | Every key, its default and its meaning, plus the two things that are easy to get wrong: `metadata_patterns`, and how the exposure levels are spaced along the fitted axis |
| [Describing your experiment](experiment-design.md) | The twelve keys that state facts about the experiment, what to write when a piece of the design is missing, the smallest configuration that runs, and what the configuration cannot say |
| [Method notes](method-notes.md) | Trimming, pooling replicates, collapsing redundant features, dropping features from the final table, and what the three figures show |
| [Reproducing the published run](reproducing.md) | The four checks against the published outputs, which numbers belong to which metadata split, where the route deliberately differs, and how to run the checks yourself |

Start with [Quickstart](quickstart.md) if you have not run the tool yet. After that:
[Configuration](configuration.md) if you are writing a run against a dataset like the published
one, and [Describing your experiment](experiment-design.md) if the design is not that one.
