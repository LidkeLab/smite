# LLM Documentation Status

Quick status report on the LLM documentation build state.

## Purpose

Provides a concise overview of:
- Build progress by phase
- Document completion status
- Pending updates
- Documentation health
- Recommendations for next actions

## Check Status File

1. Check if `.claude/llm-docs-status.json` exists
2. If missing: Report that documentation hasn't been built yet
3. If exists: Read and parse the status file
4. Present formatted status report

## Status Report Format

```
LLM Documentation Status
========================

Build Progress
--------------
Phase 1 (Essential):        ✓ Complete     (12/12 docs)   [2025-01-10]
Phase 2 (Comprehensive):    ⧗ In Progress  (18/40 docs)   [Started 2025-01-12]
Phase 3 (Maintenance):      ○ Not Started

Overall Coverage: 57.7% (30/52 planned documents)

Recent Activity
---------------
Last Build:      2025-01-10 14:30 (10 days ago)
Last Update:     2025-01-15 09:00 (5 days ago)
Last Validation: 2025-01-15 10:00 (5 days ago)

Documentation Location
----------------------
Path: doc/llm-guide/
Files: 30 documents, manifest.json, index.md

Pending Updates
---------------
⚠ 2 documents need updating:

  CRITICAL:
  • workflows/bagol-clustering.md
    Reason: Core method combineBaGoLROIs changed (API change)
    Detected: 2025-01-18 08:45 (2 days ago)
    Commits: 0a94de5, e658a5a

  MAJOR:
  • api-reference/smi-core.md
    Reason: New method added to LocalizeData
    Detected: 2025-01-19 14:20 (1 day ago)
    Commits: 999a9d2

Documentation Health
--------------------
Overall: ⚠ 95% (Good, but needs attention)

  ✓ All files exist (30/30)
  ✓ Cross-references valid (127/127)
  ✓ Frontmatter consistent (30/30)
  ✓ Internal links working (85/85)
  ⚠ Code examples not tested (15 untested)
  ⚠ 2 documents need updates

Health by Category:
  ✓ getting-started:    100% (4/4 docs current)
  ✓ core-concepts:      100% (4/4 docs current)
  ⚠ workflows:          87%  (7/8 docs current, 1 needs update)
  ✓ how-to:             100% (6/6 docs current)
  ⚠ api-reference:      75%  (3/4 docs current, 1 needs update)
  ✓ examples:           100% (3/3 docs current)
  ○ troubleshooting:    0%   (0/6 docs built)
  ○ reference:          0%   (0/5 docs built)

Recommendations
---------------
High Priority:
  1. ⚠ Apply critical updates: /smite_update_llm_docs
     2 documents need immediate attention

Medium Priority:
  2. 📝 Continue Phase 2: /smite_build_llm_docs
     22 documents remaining for complete coverage

  3. ✓ Validate updates: /smite_validate_llm_docs
     Test code examples and check consistency

Low Priority:
  4. 📊 Review documentation metrics
     Consider adding more examples to how-to guides

Access Documentation
--------------------
• Local:      file:///C:/Users/klidke/Documents/MATLAB/smite/doc/llm-guide/
• Index:      doc/llm-guide/index.md
• Manifest:   doc/llm-guide/manifest.json
• GitHub Raw: https://raw.githubusercontent.com/LidkeLab/smite/main/doc/llm-guide/
              (after pushing to repository)

Quick Commands
--------------
• Build more docs:    /smite_build_llm_docs
• Apply updates:      /smite_update_llm_docs
• Validate docs:      /smite_validate_llm_docs
• This status report: /smite_docs_status
```

## When Documentation Hasn't Been Built

If status file doesn't exist:

```
LLM Documentation Status
========================

Status: ○ Not Built

The LLM documentation has not been built yet.

Getting Started
---------------
To create interactive LLM documentation for smite:

1. Run: /smite_build_llm_docs
   This will build Phase 1 (12 essential documents)

2. Review the generated documentation
   Location: doc/llm-guide/

3. Continue with additional phases
   Run /smite_build_llm_docs again for Phase 2

What You'll Get
---------------
• Modular, LLM-friendly documentation
• Quick-start guides and tutorials
• Complete API reference
• Working code examples
• Troubleshooting guides
• Interactive user guide capabilities

Estimated Time:
• Phase 1: ~10-15 minutes (12 essential docs)
• Phase 2: ~30-45 minutes (complete coverage)
• Phase 3: Ongoing maintenance

Ready to start? Run: /smite_build_llm_docs
```

## Status Symbols

Use these symbols in the output:
- ✓ Complete / Passed / Good
- ⧗ In Progress / Running
- ○ Not Started / Planned
- ⚠ Warning / Needs Attention
- ✗ Failed / Error
- 📝 Documentation / Writing
- 📊 Metrics / Analysis
- 🔍 Validation / Testing

## Color Coding (via text descriptions)

When describing health status:
- **Excellent (95-100%):** All green, no issues
- **Good (85-94%):** Mostly green, minor warnings
- **Fair (70-84%):** Yellow, needs attention soon
- **Poor (<70%):** Red, needs immediate attention

## Detailed Statistics (Optional)

If user wants more detail, can provide:

```
Detailed Statistics
===================

Document Status Breakdown:
  Complete:      28 docs (93%)
  Needs Update:   2 docs (7%)
  In Progress:    0 docs (0%)
  Failed:         0 docs (0%)

Word Count:
  Total:         ~28,500 words
  Average:       950 words/doc
  Largest:       1,850 words (workflows/smlm-analysis.md)
  Smallest:      420 words (examples/basic-localization.md)

Last 5 Updates:
  2025-01-15 09:00 - workflows/spt-tracking.md
  2025-01-14 16:30 - api-reference/smi-namespace.md
  2025-01-14 14:15 - how-to/localize-molecules.md
  2025-01-13 11:20 - core-concepts/architecture.md
  2025-01-12 10:00 - getting-started/first-analysis.md

Git Integration:
  Derived from commits: 0a94de5, e658a5a, 999a9d2, 537a2ce, 80fa042
  Tracked source files: 45 files
  Average age: 5 days since last source change
```

## Machine-Readable Output (Optional)

If needed for automation, can output JSON:

```json
{
  "status": "in_progress",
  "coverage": 57.7,
  "health": 95.0,
  "phases": {
    "phase1": "complete",
    "phase2": "in_progress",
    "phase3": "not_started"
  },
  "documents": {
    "total": 30,
    "complete": 28,
    "needs_update": 2,
    "in_progress": 0
  },
  "pending_updates": {
    "critical": 1,
    "major": 1,
    "minor": 0
  },
  "validation": {
    "last_run": "2025-01-15T10:00:00Z",
    "files_exist": true,
    "cross_refs_valid": true,
    "code_examples_tested": false
  },
  "recommendations": [
    "Apply critical updates",
    "Continue Phase 2",
    "Validate code examples"
  ]
}
```

## Execution Notes

- Always check for status file first
- Parse dates and make them human-readable (e.g., "5 days ago")
- Sort pending updates by priority (CRITICAL > MAJOR > MINOR)
- Calculate percentages accurately
- Provide actionable recommendations
- Keep output concise but informative
- Use visual elements (symbols, spacing) for readability
- Include quick command references
- Show trends if multiple status checks (improving/degrading)

## When to Show Detailed Stats

Show detailed statistics if:
- User explicitly asks for details
- This is a milestone status check (phase complete)
- Significant changes since last check
- User wants report for documentation/records

Otherwise, keep report concise and actionable.
