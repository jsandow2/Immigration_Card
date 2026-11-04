# =============================================================================
# 00b_table_helper.R  --  write_dissertation_table()
# Part of the modular Immigration_Card_Analysis.R pipeline.
#
# EXTRACTED 2026-10-04 from 25_dissertation_table_helper.R so that 14 can use
# it too. 24 defined the function and then immediately applied it, which meant
# the house table style was only available to the script that happened to own
# it -- and 24 runs AFTER 14, so 14 could not have called it at all. The
# function now lives here and both scripts source this file. 24 keeps its
# application (the demographic summary table); nothing about the output
# changes.
#
# THE HOUSE STYLE. Matches the LaTeX table conventions used in Chapter 2
# ("Do Guns Make Us Safe?"): \caption and \label before the table, booktabs
# rules, the tabular wrapped in \resizebox{\textwidth}{!}{...}, and a manual
# \begin{minipage} footnote block below -- NOT kableExtra's built-in
# footnote()/kable_styling(), which produces a different, less polished
# structure.
#
# PREAMBLE REQUIREMENTS: booktabs, graphicx (for \resizebox), and makecell if
# any cell uses kableExtra::linebreak(). Note that this style uses NO striping,
# so it does not need \usepackage[table]{xcolor} and is not exposed to the
# xcolor/TikZ option clash.
# =============================================================================

library(kableExtra)
library(dplyr)

# data    -- data.frame/data.table to render (already formatted as strings:
#            do rounding, $-signs, comma separators, etc. before calling this)
# caption -- table caption text (LaTeX, no surrounding braces needed)
# label   -- table label, e.g. "tab:demographic_summary" (no \label{} wrapper)
# notes   -- notes text for the \textit{Notes:} block (LaTeX, plain string)
# align   -- optional column alignment string, e.g. "lccccccccc"; if NULL,
#            first column is left-aligned and the rest centered (matches
#            the guns-chapter convention of "lccccc"-style tables)
# MIDRULE SENTINEL, added 2026-10-05. Chapter 1's results tables separate the
# coefficient block from the diagnostics block (fixed effects, clustering,
# first stage, observations) with a \midrule. kable() has no facility for a
# rule in the middle of a table, so a row whose first cell is exactly
# "@MIDRULE@" is replaced by \midrule after rendering. Rows are emitted by
# kable in order, so the sentinel lands exactly where it was placed.
write_dissertation_table <- function(data, caption, label, notes,
                                     align = NULL, file = NULL) {
  if (is.null(align)) {
    align <- paste0("l", strrep("c", ncol(data) - 1))
  }

  tab_body <- kable(
    data,
    format    = "latex",
    booktabs  = TRUE,
    align     = align,
    linesep   = "",
    row.names = FALSE,
    escape    = FALSE
  )
  tab_body <- paste(
    sub("^@MIDRULE@.*$", "\\\\midrule", strsplit(as.character(tab_body), "\n")[[1]]),
    collapse = "\n")

  out <- paste0(
    "\\begin{table}[!htbp]\n",
    "\\centering\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\resizebox{\\ifdim\\width>\\textwidth \\textwidth\\else \\width\\fi}{!}{\n",
    tab_body, "\n",
    "}\n",
    "\\begin{minipage}{\\textwidth}\n",
    "\\footnotesize\n",
    "\\vspace{4pt}\n",
    "\\textit{Notes:} ", notes, "\n",
    "\\end{minipage}\n",
    "\\end{table}\n"
  )

  if (!is.null(file)) {
    writeLines(out, file)
    cat(sprintf("Table written to %s\n", file))
  }

  invisible(out)
}
