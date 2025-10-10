# Build LLM Documentation

You are tasked with building interactive LLM documentation for the smite MATLAB package.

## Overview

This command creates modular, LLM-friendly documentation in `doc/llm-guide/` that enables interactive user guidance. The documentation is built in phases and tracks progress in `.claude/llm-docs-status.json`.

## Check Current Status

1. Check if `.claude/llm-docs-status.json` exists
2. If it exists, read it to understand current build state
3. If it doesn't exist, this is a fresh build
4. Report current status to the user clearly

## Documentation Structure

Create this directory structure in `doc/llm-guide/`:

```
doc/llm-guide/
  index.md                      # Master navigation hub
  manifest.json                 # Machine-readable document catalog

  getting-started/              # Beginner entry point
    quickstart.md               # 5-minute getting started
    installation.md             # Complete setup guide
    first-analysis.md           # First complete workflow
    understanding-results.md    # Interpreting output

  core-concepts/                # Foundational knowledge
    architecture.md             # Big picture: SMF/SMD paradigm
    smf-structure.md            # Complete SMF reference
    smd-structure.md            # Complete SMD reference
    tr-structure.md             # Tracking results structure
    coordinate-system.md        # Image coordinate conventions
    data-flow.md                # Raw data → SMF → analysis → SMD

  workflows/                    # Complete analysis pipelines
    smlm-analysis.md            # Full SMLM workflow
    spt-tracking.md             # SPT workflow
    bagol-clustering.md         # BaGoL workflow
    batch-processing.md         # Publish class usage
    drift-correction.md         # Drift correction methods
    frame-connection.md         # Frame connection methods
    channel-registration.md     # Multi-channel workflows
    roi-analysis.md             # Working with ROIs

  how-to/                       # Task-based guides
    load-data.md                # Load .h5/.mat files
    create-smf.md               # Build and configure SMF
    localize-molecules.md       # LocalizeData usage
    threshold-results.md        # Apply thresholds
    visualize-results.md        # Generate SR images
    save-and-load-smd.md        # Persist results
    work-with-rois.md           # ROI extraction and analysis
    simulate-data.md            # Use simulation tools
    calibrate-camera.md         # Camera calibration
    tune-parameters.md          # Optimize SMF settings
    combine-datasets.md         # Multi-dataset analysis
    export-results.md           # Export to other formats
    use-gpu.md                  # GPU acceleration

  api-reference/                # Namespace documentation
    navigation.md               # How to read API docs
    smi-namespace.md            # Top-level classes
    smi-core.md                 # Core processing classes
    smi-cluster.md              # Clustering algorithms
    smi-stat.md                 # Statistical methods
    smi-vis.md                  # Visualization tools
    smi-sim.md                  # Simulation classes
    smi-psf.md                  # PSF modeling
    smi-helpers.md              # Utility functions

  examples/                     # Practical examples
    basic-localization.md       # Simple localization
    tracking-diffusion.md       # Track and analyze diffusion
    clustering-analysis.md      # Cluster and analyze
    multi-channel.md            # Multi-color experiment
    3d-localization.md          # Astigmatic 3D SMLM

  troubleshooting/              # Problem-solving
    installation-issues.md      # Setup problems
    gpu-problems.md             # CUDA/GPU issues
    memory-issues.md            # Out of memory errors
    compilation-errors.md       # mex/CUDA compilation
    localization-quality.md     # Poor results
    common-mistakes.md          # Frequent errors

  reference/                    # Technical details
    file-formats.md             # .h5/.mat structure
    data-structures.md          # Complete SMF/SMD/TR specs
    compilation.md              # Building mex/CUDA
    dependencies.md             # Required toolboxes
    testing.md                  # Running unit tests
```

## Build Phases

**Phase 1 - Essential (12 documents):**
Minimum viable documentation set for basic LLM interaction.

1. index.md
2. getting-started/quickstart.md
3. getting-started/installation.md
4. getting-started/first-analysis.md
5. core-concepts/architecture.md
6. core-concepts/smf-structure.md
7. core-concepts/smd-structure.md
8. workflows/smlm-analysis.md
9. workflows/spt-tracking.md
10. how-to/load-data.md
11. how-to/localize-molecules.md
12. examples/basic-localization.md

**Phase 2 - Comprehensive (30-40 documents):**
Complete all categories for full coverage.

**Phase 3 - Maintenance:**
Ongoing updates and refinements based on usage.

## Document Template

Each document MUST follow this structure with complete frontmatter:

```markdown
---
title: "Descriptive Title"
category: "getting-started|core-concepts|workflows|how-to|api-reference|examples|troubleshooting|reference"
level: "beginner|intermediate|advanced"
tags: ["tag1", "tag2", "tag3"]
prerequisites: ["doc1.md", "doc2.md"]
related: ["related1.md", "related2.md"]
summary: "One-sentence description for LLM discovery"
estimated_time: "5-10 minutes"
last_updated: "YYYY-MM-DD"
status: "complete"
---

# Title

## Purpose
What this document teaches and why it matters (1-2 paragraphs).

## Prerequisites
What you should know first. Link to prerequisite documents.

## Overview
Brief conceptual introduction (2-3 paragraphs).

## Main Content
Step-by-step instructions, explanations, code examples.
Use clear subsections.

### Code Example
\`\`\`matlab
% Working code that users can run
SMF = smi_core.SingleMoleculeFitting();
% Include comments explaining each step
\`\`\`

**Expected Output:**
```
% What users should see when they run this
```

**Explanation:**
What the code does and why.

## Common Issues
Troubleshooting specific to this topic.

## See Also
- [Related Topic 1](../path/to/doc.md) - Brief description
- [Related Topic 2](../path/to/doc.md) - Brief description
- External: [Resource Name](URL) - Brief description

## Next Steps
Suggested progression:
- [Next Topic A](../path/to/doc.md) - For going deeper
- [Next Topic B](../path/to/doc.md) - For related workflows
```

## Source Material Locations

Use these existing materials as sources:

- **Main README:** `README.md` - Overview and installation
- **Core concepts:** `doc/CoreOverview.md` - Architecture overview
- **Data structures:** `doc/DataStructures/SMF.md`, `doc/DataStructures/SMD.md`, `doc/DataStructures/TR.md`
- **Class documentation:** `MATLAB/+smi*/README.md` files in each namespace
- **Class details:** Individual `.m` files (read header comments)
- **Examples:** `MATLAB/examples/` directory
- **Extended examples:** `doc/ExtExamples/` directory
- **File formats:** `doc/FileFormats/` directory
- **Classes list:** `doc/SMITEclasses.md`

## Execution Steps

### 1. Determine Phase

- If status file doesn't exist: Start Phase 1
- If Phase 1 complete: Offer Phase 2
- If Phase 2 complete: Maintenance mode
- Ask user which phase to build (provide recommendation)

### 2. Create Directory Structure

```bash
mkdir -p doc/llm-guide/{getting-started,core-concepts,workflows,how-to,api-reference,examples,troubleshooting,reference}
```

### 3. Initialize or Update Status File

Create `.claude/llm-docs-status.json` with initial structure or update existing.

### 4. Build Each Document

For each document in the selected phase:

a. **Read source materials** relevant to the document
b. **Generate content** following the template
c. **Include working code examples** that users can run
d. **Add proper cross-references** to other docs
e. **Save the document**
f. **Update status file** with document metadata:
   - Mark document as complete
   - Record word count
   - List source files used
   - Set validation status to "passed"
   - Record current date
g. **Report progress** to user

IMPORTANT: Update status file after EACH document so work can be resumed if interrupted.

### 5. Generate manifest.json

After all documents in phase are complete:

- Parse frontmatter from all documents
- Create `doc/llm-guide/manifest.json` with:
  - Document catalog
  - Categories
  - Tags
  - Cross-references
  - Access URLs

Structure:
```json
{
  "name": "smite LLM Interactive Guide",
  "version": "1.0.0",
  "base_path": "doc/llm-guide/",
  "generated": "YYYY-MM-DDTHH:MM:SSZ",
  "index": "index.md",
  "categories": [
    {
      "id": "getting-started",
      "title": "Getting Started",
      "description": "Beginner introduction to smite",
      "document_count": 4
    }
  ],
  "documents": [
    {
      "id": "quickstart",
      "title": "Quick Start Guide",
      "path": "getting-started/quickstart.md",
      "category": "getting-started",
      "level": "beginner",
      "tags": ["setup", "first-steps"],
      "summary": "Get started with smite in 5 minutes",
      "prerequisites": [],
      "related": ["installation", "first-analysis"]
    }
  ],
  "access_methods": {
    "github_raw": "https://raw.githubusercontent.com/LidkeLab/smite/main/doc/llm-guide/",
    "local": "file:///C:/Users/klidke/Documents/MATLAB/smite/doc/llm-guide/"
  }
}
```

### 6. Generate index.md

Create the master navigation document with:
- Quick navigation paths ("I'm new", "I want to...", etc.)
- Complete document map organized by category
- Document summaries from frontmatter
- Usage instructions for both LLMs and humans

### 7. Update Status File

Mark phase as complete with:
- Completion timestamp
- Document count
- Total word count
- Coverage percentage

### 8. Final Report

Provide comprehensive report:
```
Phase 1 Build Complete!
=======================

Documents Created: 12
Total Words: ~11,000
Location: doc/llm-guide/

Files Generated:
✓ 12 documentation files
✓ manifest.json
✓ index.md
✓ Status tracking: .claude/llm-docs-status.json

Access Methods:
• Local: file:///C:/Users/klidke/Documents/MATLAB/smite/doc/llm-guide/
• GitHub Raw: https://raw.githubusercontent.com/LidkeLab/smite/main/doc/llm-guide/
  (after pushing to repository)

Next Steps:
1. Review generated documentation
2. Test by asking LLM questions about smite
3. Run /smite_validate_llm_docs to check quality
4. Continue with Phase 2: /smite_build_llm_docs
5. Or check status anytime: /smite_docs_status

Documentation Coverage: 23% (12/52 planned documents)
```

## Content Guidelines

**Practical:**
- Every tutorial must have working code examples
- Code should be complete and runnable
- Show expected output
- Explain what the code does

**Progressive:**
- Start with simple examples
- Build complexity gradually
- Link to prerequisites clearly
- Suggest next steps

**Connected:**
- Link related documents
- Reference source code files with line numbers when relevant
- Cross-reference between sections
- Build knowledge progressively

**Complete:**
- Answer: What is this?
- Answer: Why does it matter?
- Answer: How do I use it?
- Answer: When should I use it?
- Include troubleshooting

**Accurate:**
- Test code examples
- Verify against current source code
- Check parameter names and values
- Ensure paths and references are correct

## Important Notes

- Create `.claude/llm-docs-status.json` if it doesn't exist
- Update status after EACH document for resumability
- If interrupted, can resume from status file
- Keep manifest.json in sync with actual documents
- Don't regenerate existing docs unless user requests
- For Phase 2+, consult status file for what's already done
- Always provide clear progress reporting
- Each document should be self-contained but connected

## Error Handling

If you encounter issues:
- Report what failed and why
- Save progress in status file before stopping
- Suggest how to fix the issue
- Allow user to resume after fixes

## Validation

Before marking phase complete:
- Verify all documents exist
- Check all internal cross-references resolve
- Ensure manifest matches reality
- Confirm frontmatter is complete
- Save validation results to status file
