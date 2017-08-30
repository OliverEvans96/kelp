# Kelp Project Phase 1 Outline

---

### Goal
Given kelp configuration & light state, determine how
much light is absorbed by kelp at each depth.


### Inputs:
- ind(z)
- A(z)
- IOPs
- Light Surface BC

### Outputs:
- lightAbsorbed(z)

### Steps:

1. Generate kelp
- RTE to calculate light field
- Integrate light over kelp to calculate absorption