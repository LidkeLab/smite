### +smi_core/@TrackingResults

TrackingResults: A class defining the Tracking Results structure

This datatype is one of the primary results structures in the smite
enviroment. The Tracking Results (TR) structure is an input/output of
many methods in smite which are related to single-particle tracking
(SPT) analysis.  TR structures are organized as follows: each
structure element corresponds to a single trajectory, i.e., TR(n)
contains all relevant properties of the n-th trajectory. The TR
structure is intended to carry all necessary information from an SMD
structure (see (smi_core.SingleMoleculeData)[SMD.md]) but organized
in a more user-friendly manner for SPT data.

The TR structure is just an array of SMD structures, with each array
element being an SMD structure corresponding to a single trajectory.

SEE ALSO:
- [smi_core.SingleMoleculeData](../@SingleMoleculeData/README.md)
- [smi_core.SingleMoleculeFitting](../@SingleMoleculeFitting/README.md)

---

methods:
- **[catTR](catTR.m)**:
  concatenates two TR structures
- **[computeTrajDurations](computeTrajDurations.m)**:
  computes the duration of trajectories in TR
- **[computeTrajFidelity](computeTrajFidelity.m)**:
  computes the trajectory length divided by duration
- **[computeTrajLengths](computeTrajLengths.m)**:
  computes the length of trajectories in TR
- **[convertSMDToTR](convertSMDToTR.m)**:
  converts an SMD into a TR structure
- **[convertTRToSMD](convertTRToSMD.m)**:
  converts a TR back into an SMD structure.
- **[getTRIndex](getTRIndex.m)**:
  returns the indices in TR corresponding to TrajectoryIDs
- **[joinTraj](joinTraj.m)**:
  joins a set of trajectories into one trajectory
- **[padTR](padTR.m)**:
  pads the input 'TR' with fields present in 'TRPadding'
- **[threshTrajLength](threshTrajLength.m)**:
  removes trajectories smaller than a minimum track length
