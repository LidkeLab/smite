---
title: "Navigating the API Reference"
category: "api-reference"
level: "beginner"
tags: ["api", "documentation", "reference", "methods", "properties", "reading-docs"]
prerequisites: ["../getting-started/quickstart.md", "../core-concepts/architecture.md"]
related: ["smi-namespace.md", "smi-core.md", "../how-to/create-smf.md"]
summary: "Learn how to read and navigate smite's API reference documentation effectively, including understanding method signatures, property tables, and cross-references"
estimated_time: "15 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# Navigating the API Reference

## Purpose

The API reference documentation provides complete technical specifications for every class, method, and property in smite. This guide teaches you how to read and navigate these references efficiently, helping you find exactly what you need and understand technical specifications without getting overwhelmed.

## Prerequisites

- Basic MATLAB experience
- Familiarity with [smite architecture](../core-concepts/architecture.md)
- Understanding of [basic workflow concepts](../getting-started/quickstart.md)

## Overview

Reading API documentation is a skill that improves with practice. This guide covers:

1. **Organization** - How API references are structured
2. **Reading Method Signatures** - Decoding function call syntax
3. **Properties vs Methods** - Understanding the distinction
4. **Static vs Instance Methods** - When to use each
5. **Parameter Tables** - Interpreting inputs and outputs
6. **Finding Examples** - Locating practical usage code
7. **Cross-Referencing** - Navigating between related documentation

## API Reference Organization

### Namespace Structure

smite's API documentation mirrors the code organization using MATLAB's `+` namespace system:

```
+smi              Top-level workflows (SMLM, SPT, BaGoL, Publish)
+smi_core         Core processing (LocalizeData, DriftCorrection, etc.)
+smi_sim          Simulation tools
+smi_cluster      Clustering algorithms
+smi_stat         Statistical methods
+smi_vis          Visualization tools
+smi_psf          Point spread function tools
+smi_helpers      Helper utilities
```

**Navigation strategy:**

- **Start with +smi** for complete workflows and main entry points
- **Use +smi_core** when you need fine-grained control over processing steps
- **Explore specialized namespaces** (+smi_cluster, +smi_stat) for specific analyses
- **Refer to +smi_vis** when generating figures and visualizations

### Document Layout

Each API reference document follows a consistent structure:

```markdown
# API Reference: +namespace

## Overview
[High-level description of the namespace]

## Class Reference
### ClassName1
  - Description
  - Constructor
  - Properties
  - Methods
  - Usage Examples

### ClassName2
  ...
```

**Reading tip:** Scan the Overview section first to understand what the namespace provides, then jump directly to the class you need.

## Understanding Method Signatures

Method signatures show you exactly how to call a function. Let's break down the notation.

### Basic Signature Format

```matlab
[Output1, Output2] = obj.methodName(Input1, Input2, Input3)
```

**Components:**

- `[Output1, Output2]` - What the method returns (left side of `=`)
- `obj` - The object instance you're calling the method on
- `.methodName` - The method name
- `(Input1, Input2, Input3)` - Required inputs inside parentheses

### Example: Reading a Real Signature

From the LocalizeData class:

```matlab
[SMD, SMDPreThresh] = LD.genLocalizations()
```

**Interpretation:**

- **LD** - Your LocalizeData object instance
- **genLocalizations()** - The method name (no inputs required)
- **[SMD, SMDPreThresh]** - Returns two outputs:
  - `SMD` - Thresholded localizations
  - `SMDPreThresh` - Unthresholded localizations

**Usage:**

```matlab
LD = smi_core.LocalizeData(Data, SMF);
[SMD, SMDPreThresh] = LD.genLocalizations();
```

### Optional Parameters

When inputs are optional, documentation shows multiple signature variations:

```matlab
obj = smi.SMLM()              % Opens GUI
obj = smi.SMLM(SMF)           % Uses provided SMF, no GUI
obj = smi.SMLM(SMF, StartGUI) % Control GUI opening
```

**Interpretation:**

- First line: No inputs, uses all defaults
- Second line: One input (SMF), other parameters use defaults
- Third line: Two inputs, full control

**Reading tip:** Try the simplest signature first (fewest inputs), then add parameters as needed.

### Variable Arguments

Some methods accept variable numbers of inputs:

```matlab
[~, Data, SMF] = LD.loadRawData(SMF, varargin)
```

The `varargin` notation means "additional optional arguments." Check the method's parameter table to see what these can be.

## Properties vs Methods

Understanding the difference between properties and methods is essential for using classes correctly.

### Properties (Data Storage)

Properties hold data about an object. They're like variables attached to the object.

**Access syntax:**

```matlab
value = obj.PropertyName         % Get property value
obj.PropertyName = newValue      % Set property value
```

**Example from SMLM class:**

```matlab
SMLMobj.Verbose = 2;             % Set verbosity level
zoom = SMLMobj.SRImageZoom;      % Get current zoom value
```

**Properties are nouns** - they describe characteristics (Verbose, SRImageZoom, SMD, SMF).

### Methods (Actions)

Methods perform operations. They're like functions attached to an object.

**Call syntax:**

```matlab
output = obj.methodName(inputs)  % Call method with inputs
obj.methodName()                 % Call method with no inputs
```

**Example from SMLM class:**

```matlab
SMLMobj.fullAnalysis();          % Perform complete analysis
SMLMobj.testFit(1);              % Test on dataset 1
SMLMobj.generatePlots(dir1, dir2, 'test', true, plots);
```

**Methods are verbs** - they do things (fullAnalysis, testFit, generatePlots).

### Property vs Method: How to Tell

In API documentation:

**Properties appear in tables:**

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `Verbose` | integer | 1 | Verbosity level |
| `SMD` | struct | [] | Results structure |

**Methods appear with signatures:**

```matlab
obj.fullAnalysis()
obj.testFit(DatasetIndex)
```

**Quick test in MATLAB:**

```matlab
% Properties don't use parentheses
obj.Verbose       % Property - no ()

% Methods require parentheses (even if empty)
obj.fullAnalysis()  % Method - needs ()
```

## Static vs Instance Methods

Methods come in two flavors: static and instance.

### Instance Methods

Instance methods operate on a specific object. You need to create an object first.

**Pattern:**

```matlab
obj = ClassName();          % Create instance
result = obj.methodName();  % Call method on instance
```

**Example:**

```matlab
% Create SMLM object
SMLMobj = smi.SMLM(SMF);

% Call instance method
SMLMobj.fullAnalysis();     % Operates on this specific SMLMobj
```

**Instance methods often modify the object** or use its properties.

### Static Methods

Static methods don't require an object instance. They're utility functions attached to a class.

**Pattern:**

```matlab
result = ClassName.methodName(inputs);  % Call directly on class
```

**Example:**

```matlab
% No object creation needed
counts = smi.SMLM.fitsPerFrame(SMD);

% Another example
SMDCombined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2);
```

**Static methods are pure functions** - they take inputs and return outputs without modifying any object state.

### How to Identify Static Methods

In API documentation, static methods are labeled:

```markdown
### Static Methods

##### methodName
[Static method description]
```

Or shown with class name instead of `obj`:

```matlab
% Instance method signature
[SMD, SMDPreThresh] = obj.genLocalizations()

% Static method signature
SMD = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2)
```

**Reading tip:** If you see `ClassName.method` instead of `obj.method`, it's static.

## Reading Parameter Tables

Parameter tables specify what inputs a method expects and what outputs it returns.

### Input Parameter Tables

**Format:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `SMF` | struct | required | Single Molecule Fitting structure |
| `StartGUI` | logical | true | Open GUI flag |
| `Verbose` | integer | 1 | Verbosity level |

**Interpretation:**

- **Parameter column** - Name you use in code
- **Type column** - What kind of data (struct, logical, integer, array, etc.)
- **Default column** - What happens if you don't provide it
  - "required" means you must provide it
  - A value means it's optional with that default
- **Description column** - What the parameter does

**Example usage:**

```matlab
% Using all defaults where possible
obj = smi.SPT(SMF);              % Only required parameter

% Overriding defaults
obj = smi.SPT(SMF, false);       % Set StartGUI = false
```

### Output Tables

Some methods return multiple values:

**Returns:**

| Output | Type | Description |
|--------|------|-------------|
| `TR` | struct array | Tracking Results |
| `SMD` | struct | Single Molecule Data |
| `SMDPreThresh` | struct | Pre-threshold SMD |

**Capturing outputs:**

```matlab
% Capture all outputs
[TR, SMD, SMDPreThresh] = SPTobj.performFullAnalysis();

% Capture only what you need
[TR, SMD] = SPTobj.performFullAnalysis();  % Skip SMDPreThresh

% Capture just the first
TR = SPTobj.performFullAnalysis();

% Ignore specific outputs with ~
[TR, ~, SMDPreThresh] = SPTobj.performFullAnalysis();  % Skip SMD
```

**Reading tip:** You can always capture fewer outputs than documented, but MATLAB returns them in the order listed.

## Finding Usage Examples

API documentation includes examples at multiple locations.

### Where Examples Appear

1. **After each method description:**

```matlab
##### methodName

Description of what the method does.

**Example:**
```matlab
obj = ClassName();
result = obj.methodName(input);
```
```

2. **In "Usage Patterns" sections:**

Near the end of each class reference, find realistic multi-step workflows:

```markdown
### Usage Patterns

#### Pattern 1: Basic Usage
[Complete working example]

#### Pattern 2: Advanced Usage
[More complex scenario]
```

3. **In "See Also" sections:**

Links to related how-to guides with extended examples:

```markdown
### See Also
- [How to Localize Molecules](../how-to/localize-molecules.md)
- [SMLM Workflow](../workflows/smlm-analysis.md)
```

### How to Use Examples Effectively

**Start with the simplest example:**

```matlab
% From API reference:
LD = smi_core.LocalizeData(Data, SMF);
[SMD, ~] = LD.genLocalizations();
```

**Modify for your needs:**

```matlab
% Your version with your data
LD = smi_core.LocalizeData(myData, mySMF);
LD.Verbose = 2;  % Add more output
[SMD, SMDPreThresh] = LD.genLocalizations();  % Keep both outputs
```

**Combine multiple examples:**

API references show individual methods. Combine them into workflows:

```matlab
% Combine LoadData + DataToPhotons + LocalizeData examples
LD = smi_core.LoadData();
[~, RawData, SMF] = LD.loadRawData(SMF, 1);

DTP = smi_core.DataToPhotons(SMF, RawData);
[Data, ~] = DTP.convertData();

LD = smi_core.LocalizeData(Data, SMF);
[SMD, ~] = LD.genLocalizations();
```

## Cross-Referencing Between Documents

API documentation is interconnected. Learning to navigate these connections saves time.

### Types of Cross-References

**1. Prerequisites (at top of document):**

```yaml
prerequisites: ["../core-concepts/architecture.md", "../core-concepts/smf-structure.md"]
```

Read these first if you're confused by concepts in the current document.

**2. Related documents:**

```yaml
related: ["../workflows/smlm-analysis.md", "../how-to/localize-molecules.md"]
```

Read these for complementary information (workflows show how to use APIs together).

**3. In-text "See Also" sections:**

```markdown
**See Also:**
- [SMD Structure Guide](../core-concepts/smd-structure.md)
- [How to Threshold Results](../how-to/threshold-results.md)
```

Follow these for detailed explanations of related topics.

### Navigation Strategy

**Top-down approach (recommended for beginners):**

1. Start with [Architecture Overview](../core-concepts/architecture.md)
2. Read [+smi namespace](smi-namespace.md) for main workflows
3. Dive into [+smi_core namespace](smi-core.md) for implementation details
4. Explore specialized namespaces as needed

**Bottom-up approach (when you know what you need):**

1. Use your editor's search to find method names
2. Read that specific method's documentation
3. Check "See Also" links if you need more context
4. Refer to prerequisites only if you're confused

**Problem-solving approach:**

1. Start with [How-to guides](../how-to/) for task-based help
2. Follow API reference links from how-to guides
3. Use examples from API reference
4. Return to how-to guide to continue your task

## Practical Reading Examples

Let's practice reading real API documentation.

### Example 1: Using LocalizeData

**Goal:** Understand how to use the LocalizeData class.

**Step 1: Find the class**

Navigate to: [+smi_core namespace](smi-core.md) → LocalizeData section

**Step 2: Read the description**

> LocalizeData performs complete localization pipeline from photon-converted images to thresholded localizations.

**Translation:** This class takes your camera data (after converting to photons) and finds molecules.

**Step 3: Check the constructor**

```matlab
LD = smi_core.LocalizeData(ScaledData, SMF, Verbose, ResultsDir)
```

**Translation:** You need at minimum `ScaledData` and `SMF`. The others are optional.

**Step 4: Find the main method**

```matlab
[SMD, SMDPreThresh] = LD.genLocalizations()
```

**Translation:** Call this method to get results. It returns two versions: filtered (SMD) and unfiltered (SMDPreThresh).

**Step 5: Try the example**

```matlab
% From documentation
LD = smi_core.LocalizeData(Data, SMF);
[SMD, SMDPreThresh] = LD.genLocalizations();

% Adapt to your needs
LD = smi_core.LocalizeData(myPhotonData, mySMF, 2);  % Verbose=2
[results, allFits] = LD.genLocalizations();
fprintf('Found %d molecules\n', length(results.X));
```

### Example 2: Understanding Property Tables

**From smi.SMLM documentation:**

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `Verbose` | integer | 1 | Verbosity level (0=silent, 1=normal, 2+=debug) |
| `SRImageZoom` | numeric | 20 | Magnification factor for super-resolution images |

**Reading this table:**

```matlab
% Access properties
SMLMobj = smi.SMLM(SMF);

% Get current value
currentVerbose = SMLMobj.Verbose;  % Returns 1 (the default)

% Set new value
SMLMobj.Verbose = 2;       % Enable debug output
SMLMobj.SRImageZoom = 50;  % Higher magnification

% Type information tells you what's valid
SMLMobj.Verbose = 0;       % OK - integer
SMLMobj.Verbose = 3;       % OK - integer
SMLMobj.Verbose = 'high';  % ERROR - must be integer, not string
```

### Example 3: Static vs Instance Methods

**From smi.SPT documentation:**

**Instance method:**

```matlab
[TR, SMD] = obj.performFullAnalysis()
```

**Static method:**

```matlab
SMD = smi.SPT.genTrajFF(SMD, SMF, RhoOff, NonLinkMarker)
```

**Using each:**

```matlab
% Instance method - need object first
SPTobj = smi.SPT(SMF);
[TR, SMD] = SPTobj.performFullAnalysis();  % obj.method()

% Static method - no object needed
SMD_linked = smi.SPT.genTrajFF(SMD, SMF, 0.1, -1);  % Class.method()
```

**When to use each:**

- **Instance methods** for complete workflows (performFullAnalysis)
- **Static methods** for utility functions you might call standalone (genTrajFF)

## Common Reading Mistakes and Solutions

### Mistake 1: Forgetting Parentheses on Methods

**Wrong:**

```matlab
SMLMobj.fullAnalysis  % Missing ()
```

**Error message:** "Reference to non-existent field 'fullAnalysis'"

**Correct:**

```matlab
SMLMobj.fullAnalysis()  % Methods always need ()
```

**Rule:** If it's a method, use `()` even if there are no inputs.

### Mistake 2: Using Instance Method Statically

**Wrong:**

```matlab
smi.SMLM.fullAnalysis()  % fullAnalysis is instance method
```

**Error message:** "Not enough input arguments"

**Correct:**

```matlab
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis()  % Need object instance
```

**Rule:** Check if method is labeled "Static" in documentation.

### Mistake 3: Wrong Number of Outputs

**Wrong:**

```matlab
result = LD.genLocalizations()  % Returns 2 outputs, capturing only 1
```

**Issue:** Not an error, but you're losing information. `result` contains only SMD, not SMDPreThresh.

**Better:**

```matlab
[SMD, SMDPreThresh] = LD.genLocalizations()  % Capture both
% or
[SMD, ~] = LD.genLocalizations()  % Explicitly ignore second output
```

**Rule:** Check "Returns" section to see all available outputs.

### Mistake 4: Required vs Optional Parameters

**Wrong:**

```matlab
obj = smi.SPT()  % Missing required SMF parameter
```

**Error message:** Varies depending on class implementation.

**Correct:**

```matlab
SMF = smi_core.SingleMoleculeFitting();
obj = smi.SPT(SMF)  % Provide required parameter
```

**Rule:** Parameter table shows "required" vs default values. Required parameters come first in signature.

## Quick Reference Card

Keep this nearby when reading API docs:

### Method Signatures

```matlab
output = obj.method(input)        % Instance method
output = Class.method(input)      % Static method
[out1, out2] = obj.method()       % Multiple outputs
obj.method()                      % No outputs (modifies object)
```

### Property Access

```matlab
value = obj.Property              % Get property
obj.Property = newValue           % Set property
```

### Reading Parameter Tables

| Symbol | Meaning |
|--------|---------|
| required | You must provide this parameter |
| value | Optional, uses this default if omitted |
| [] | Optional, empty by default |
| - | Required (no default available) |

### Documentation Sections

| Section | Contains | When to Read |
|---------|----------|--------------|
| Overview | High-level purpose | Always start here |
| Constructor | How to create objects | When instantiating classes |
| Properties | Object data fields | When configuring objects |
| Methods | What operations are available | When using the class |
| Usage Patterns | Complete workflows | When learning the class |
| See Also | Related documentation | When you need more context |

## Next Steps

Now that you understand how to read API documentation:

1. **Practice:** Pick a class you're interested in (e.g., [smi.SMLM](smi-namespace.md#smlm)) and read through it
2. **Try examples:** Run the code examples from documentation in MATLAB
3. **Explore:** Follow "See Also" links to build a mental map of how classes relate
4. **Apply:** Use API references while following [how-to guides](../how-to/) for specific tasks

## See Also

### API References
- [+smi Namespace](smi-namespace.md) - Main workflow classes
- [+smi_core Namespace](smi-core.md) - Core processing classes
- [+smi_vis Namespace](smi-vis.md) - Visualization tools

### Conceptual Guides
- [Architecture Overview](../core-concepts/architecture.md) - How smite is organized
- [SMF Structure](../core-concepts/smf-structure.md) - Parameter structure reference
- [SMD Structure](../core-concepts/smd-structure.md) - Results structure reference

### Practical Guides
- [How to Create SMF](../how-to/create-smf.md) - Configuring analysis parameters
- [How to Localize Molecules](../how-to/localize-molecules.md) - Using core classes
- [SMLM Workflow](../workflows/smlm-analysis.md) - Putting it all together

---

**Reading API documentation is a skill.** Don't worry if it feels overwhelming at first. Start with high-level workflow classes (smi.SMLM, smi.SPT), use the examples, and gradually explore more detailed references as you need them. The documentation is comprehensive so you can find answers—but you don't need to read everything at once.
