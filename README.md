# Hurdle Structure Latent Class Mixed Model (HS-LCMM)
This is the base code of the article "Modeling New Product Adoption: A Hurdle Structure Latent Class Mixed Model Applied to Plant-Based Meat Alternatives". We follow: "https://github.com/CecileProust-Lima/lcmm" for the Second Hurdle part. It is a very flexible model and has room for customization.
   
1. **Second‑hurdle estimation**  
   The implementation follows the methodology of the **LCMM** R package, using its `hlme()` function by default for the second‑hurdle mixed model.

2. **Customizing the joint log‑likelihood**  
   Our reference code maximizes a joint likelihood that includes *only fixed effects,* mirroring the setup in the accompanying paper.  
   If you wish to incorporate random effects (or any other extensions), you will need to modify the LCMM log‑likelihood accordingly.  
   Detailed guidance and examples can be found in the LCMM package documentation and vignette.

3. **Semi‑parametric alternative**  
   Users who prefer a semi‑parametric specification can simply switch from `hlme()` to `lcmm()`. All other pieces of the workflow remain unchanged.

4. **Tip for large time variables**  
   When the time index is a large integer (e.g., measured in seconds, minutes, or age), nonlinear terms in `lcmm()` may become numerically unstable. In practice, rescaling the time variable—dividing by 10 to 100—improves convergence and guards against sudden crashes without affecting the substantive interpretation of coefficients. 
