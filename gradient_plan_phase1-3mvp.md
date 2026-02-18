# Plan: Add Graded Junction Support to P-N Junction Simulator

## Context

The current P-N junction simulator models **abrupt junctions** with constant doping on each side (N_D on n-side, N_A on p-side) and a sharp transition at x=0. This limits its educational and research value, as real semiconductor junctions often have **graded doping profiles** from diffusion, ion implantation, or epitaxial growth.

This plan adds support for **graded junctions** with spatially-varying doping N_D(x) and N_A(x), enabling users to explore how different doping profiles affect:
- Depletion width and electric field distribution
- Built-in potential and band bending
- Junction capacitance and breakdown voltage

The implementation uses a **full numerical Poisson solver** to accurately model arbitrary doping profiles while maintaining backward compatibility with the existing abrupt junction mode.

## User Requirements

Based on user responses:
1. **Profile types:** Linear, Exponential, Gaussian, and Custom/Arbitrary expressions
2. **Physics solver:** Full numerical Poisson solver (most accurate)
3. **UI approach:** Advanced mode toggle (preserves simple mode for beginners)

## Critical File

**`/Users/vblake/isd/pneg/pn-junction-simulator.html`** - Single-file application (all changes in this file)

## Implementation Overview

### Phase 1: Core Numerical Solver (MVP)
Build the foundation for graded junctions with linear profiles only.

### Phase 2: UI Integration
Connect solver to existing visualization charts.

### Phase 3: Band Diagram Integration
Use numerical φ(x) for accurate band bending.

### Phase 4: Additional Profiles
Add exponential, gaussian, and custom expressions.

### Phase 5: Visualization Enhancements
Add doping profile chart and improved info displays.

### Phase 6: Polish & Testing
Error handling, validation, edge cases, optimization.

---

## Detailed Implementation Steps

## PHASE 1: Core Numerical Solver (Lines 445-540)

### Step 1.1: Add Data Structures
**Location:** After line 444 (after `currentMaterial` declaration)

Add configuration objects:
```javascript
// Graded junction configuration
let junctionConfig = {
    mode: 'abrupt',  // 'abrupt' or 'graded'
    profileType: 'linear',  // 'linear', 'exponential', 'gaussian', 'custom'
    nProfile: {
        type: 'constant',
        N0: null,  // Set from nConc slider
        params: { a: 1e20 }  // Linear gradient (cm^-4)
    },
    pProfile: {
        type: 'constant',
        N0: null,
        params: { a: -1e20 }
    }
};

// Grid configuration
const GRID_CONFIG = {
    points: 500,
    depletionFactor: 0.3  // Extra margin beyond depletion
};

// Grid data storage
let gridData = {
    x: [],      // Position (m)
    ND: [],     // N_D(x) (m^-3)
    NA: [],     // N_A(x) (m^-3)
    rho: [],    // Charge density (C/m^3)
    phi: [],    // Potential (V)
    E: [],      // E-field (V/m)
    xp: 0,      // P-side depletion edge
    xn: 0       // N-side depletion edge
};
```

### Step 1.2: Implement Profile Evaluation
**Location:** After `calculateDepletionParams()` (after line 536)

```javascript
function evaluateProfile(profile, x_cm) {
    // x_cm in centimeters, returns concentration in cm^-3
    if (profile.type === 'constant') {
        return profile.N0;
    }

    const params = profile.params;

    switch (profile.type) {
        case 'linear':
            // N(x) = N0 + a*x
            return Math.max(0, profile.N0 + params.a * x_cm);

        case 'exponential':
            // N(x) = N0 * exp(lambda*x)
            return profile.N0 * Math.exp(params.lambda * x_cm);

        case 'gaussian':
            // N(x) = N0 * exp(-(x/sigma)^2)
            const arg = x_cm / params.sigma;
            return profile.N0 * Math.exp(-arg * arg);

        default:
            return profile.N0;
    }
}
```

### Step 1.3: Implement Grid Creation
```javascript
function createAdaptiveGrid(ND_profile, NA_profile, Vbi, eps) {
    // Estimate depletion width
    const ND_avg = evaluateProfile(ND_profile, 0);
    const NA_avg = evaluateProfile(NA_profile, 0);
    const W_est = Math.sqrt(2 * eps * Vbi / q * (1/(NA_avg*1e6) + 1/(ND_avg*1e6)));

    const xp_est = W_est * ND_avg / (NA_avg + ND_avg);
    const xn_est = W_est * NA_avg / (NA_avg + ND_avg);
    const margin = GRID_CONFIG.depletionFactor;

    const xMin = -xp_est * (1 + margin);
    const xMax = xn_est * (1 + margin);
    const N = GRID_CONFIG.points;

    const x = [];
    const dx = (xMax - xMin) / (N - 1);
    for (let i = 0; i < N; i++) {
        x.push(xMin + i * dx);
    }

    return {
        x: x,
        ND: new Array(N).fill(0),
        NA: new Array(N).fill(0),
        rho: new Array(N).fill(0),
        phi: new Array(N).fill(0),
        E: new Array(N).fill(0),
        xp: xp_est,
        xn: xn_est
    };
}
```

### Step 1.4: Implement Poisson Solver (FDM)
```javascript
function solvePoissonFDM(x, rho, eps, Vbi_eff) {
    // Solve d²φ/dx² = -ρ/ε using Finite Difference Method
    const N = x.length;

    // Tridiagonal matrix coefficients
    const a = new Array(N).fill(0);  // sub-diagonal
    const b = new Array(N).fill(0);  // diagonal
    const c = new Array(N).fill(0);  // super-diagonal
    const d = new Array(N).fill(0);  // RHS

    // Boundary conditions
    b[0] = 1; c[0] = 0; d[0] = 0;           // φ(x_min) = 0
    b[N-1] = 1; a[N-1] = 0; d[N-1] = Vbi_eff;  // φ(x_max) = Vbi

    // Interior points
    for (let i = 1; i < N - 1; i++) {
        const dx_left = x[i] - x[i-1];
        const dx_right = x[i+1] - x[i];
        const dx_avg = (dx_left + dx_right) / 2;

        a[i] = 1 / (dx_left * dx_avg);
        b[i] = -1/(dx_left*dx_avg) - 1/(dx_right*dx_avg);
        c[i] = 1 / (dx_right * dx_avg);
        d[i] = -rho[i] / eps;
    }

    return solveTridiagonal(a, b, c, d);
}

function solveTridiagonal(a, b, c, d) {
    // Thomas algorithm for tridiagonal systems
    const N = b.length;
    const c_star = new Array(N);
    const d_star = new Array(N);
    const x = new Array(N);

    // Forward sweep
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (let i = 1; i < N; i++) {
        const denom = b[i] - a[i] * c_star[i-1];
        c_star[i] = c[i] / denom;
        d_star[i] = (d[i] - a[i] * d_star[i-1]) / denom;
    }

    // Back substitution
    x[N-1] = d_star[N-1];
    for (let i = N - 2; i >= 0; i--) {
        x[i] = d_star[i] - c_star[i] * x[i+1];
    }

    return x;
}
```

### Step 1.5: Implement E-field Calculation
```javascript
function computeElectricField(x, phi) {
    // E = -dφ/dx using central difference
    const N = x.length;
    const E = new Array(N);

    // Interior points (central difference)
    for (let i = 1; i < N - 1; i++) {
        E[i] = -(phi[i+1] - phi[i-1]) / (x[i+1] - x[i-1]);
    }

    // Boundaries (forward/backward difference)
    E[0] = -(phi[1] - phi[0]) / (x[1] - x[0]);
    E[N-1] = -(phi[N-1] - phi[N-2]) / (x[N-1] - x[N-2]);

    return E;
}
```

### Step 1.6: Implement Main Solver Function
```javascript
function solvePoisson(calc) {
    const ND_profile = junctionConfig.nProfile;
    const NA_profile = junctionConfig.pProfile;
    const eps = eps0 * eps_r[currentMaterial];
    const Vbi_eff = Math.max(calc.Vbi - calc.forwardBias, 0.001);

    try {
        // Create grid
        const grid = createAdaptiveGrid(ND_profile, NA_profile, Vbi_eff, eps);

        // Evaluate doping profiles
        for (let i = 0; i < grid.x.length; i++) {
            const x_cm = grid.x[i] * 100;  // m → cm
            grid.ND[i] = evaluateProfile(ND_profile, x_cm) * 1e6;  // cm^-3 → m^-3
            grid.NA[i] = evaluateProfile(NA_profile, x_cm) * 1e6;
            grid.rho[i] = q * (grid.ND[i] - grid.NA[i]);
        }

        // Solve Poisson equation
        grid.phi = solvePoissonFDM(grid.x, grid.rho, eps, Vbi_eff);

        // Compute E-field
        grid.E = computeElectricField(grid.x, grid.phi);

        // Find depletion edges (where |ρ| drops to ~1% of peak)
        const junctionIdx = grid.x.findIndex(x => x >= 0);
        const rhoMax = Math.max(...grid.rho.map(Math.abs));
        const threshold = rhoMax * 0.01;

        for (let i = 0; i < junctionIdx; i++) {
            if (Math.abs(grid.rho[i]) < threshold) {
                grid.xp = Math.abs(grid.x[i]);
                break;
            }
        }

        for (let i = grid.x.length - 1; i > junctionIdx; i--) {
            if (Math.abs(grid.rho[i]) < threshold) {
                grid.xn = grid.x[i];
                break;
            }
        }

        return grid;

    } catch (error) {
        console.error('Poisson solver failed:', error);
        alert('Numerical solver error. Reverting to abrupt junction mode.');
        document.getElementById('advancedMode').checked = false;
        junctionConfig.mode = 'abrupt';
        return null;
    }
}
```

---

## PHASE 2: UI Integration (Lines 300-325, 1350-1380)

### Step 2.1: Add Advanced Mode Checkbox
**Location:** After line 302 (after "Show Intrinsic Level" checkbox)

```html
<div class="control-group checkbox-group">
    <input type="checkbox" id="advancedMode">
    <label for="advancedMode">⚙️ Advanced: Graded Junctions</label>
</div>
```

### Step 2.2: Add Graded Controls Section
**Location:** After main controls div (after line 321)

```html
<div id="gradedControls" class="controls" style="display: none; background: #fef3c7; border: 2px solid #f59e0b;">
    <h4 style="grid-column: 1 / -1; color: #92400e;">
        📊 Graded Junction Configuration
    </h4>

    <div class="control-group">
        <label for="profileType">Profile Type</label>
        <select id="profileType">
            <option value="linear">Linear: N(x) = N₀ + a·x</option>
            <option value="exponential">Exponential: N(x) = N₀·e^(λx)</option>
            <option value="gaussian">Gaussian: N(x) = N₀·e^(-(x/σ)²)</option>
            <option value="custom">Custom Expression</option>
        </select>
    </div>

    <!-- Linear parameters -->
    <div id="linearParamsN" class="control-group">
        <label>N-side gradient (a): <span id="linearNValue">1.00e+20</span> cm⁻⁴</label>
        <input type="range" id="linearN" min="18" max="22" step="0.1" value="20">
    </div>

    <div id="linearParamsP" class="control-group">
        <label>P-side gradient (a): <span id="linearPValue">-1.00e+20</span> cm⁻⁴</label>
        <input type="range" id="linearP" min="-22" max="-18" step="0.1" value="-20">
    </div>

    <!-- Parameters for other profiles (exponential, gaussian, custom) -->
    <!-- Similar structure, initially hidden -->

    <div class="control-group" style="grid-column: 1 / -1;">
        <div style="background: #dbeafe; padding: 10px; border-radius: 4px; font-size: 0.875rem;">
            <strong>Note:</strong> N_D and N_A sliders set base concentration N₀.
            Profile parameters define spatial variation.
        </div>
    </div>
</div>
```

### Step 2.3: Add Event Listeners
**Location:** After existing event listeners (after line 1330)

```javascript
// Advanced mode toggle
document.getElementById('advancedMode').addEventListener('change', function() {
    junctionConfig.mode = this.checked ? 'graded' : 'abrupt';
    document.getElementById('gradedControls').style.display =
        this.checked ? 'grid' : 'none';
    update();
});

// Profile type selector
document.getElementById('profileType').addEventListener('change', function() {
    junctionConfig.profileType = this.value;

    // Show/hide appropriate parameter controls
    ['linearParamsN', 'linearParamsP'].forEach(id => {
        document.getElementById(id).style.display =
            this.value === 'linear' ? 'flex' : 'none';
    });
    // Similar for exponential, gaussian, custom

    update();
});

// Parameter sliders
document.getElementById('linearN').addEventListener('input', function() {
    const value = Math.pow(10, parseFloat(this.value));
    document.getElementById('linearNValue').textContent = value.toExponential(2);
    junctionConfig.nProfile.params.a = value;
    update();
});

document.getElementById('linearP').addEventListener('input', function() {
    const value = Math.pow(10, parseFloat(this.value));
    document.getElementById('linearPValue').textContent = value.toExponential(2);
    junctionConfig.pProfile.params.a = value;
    update();
});
```

### Step 2.4: Modify calculate() Function
**Location:** Line 483, update return statement

```javascript
function calculate() {
    // ... existing code ...

    // NEW: Solve graded junction if in advanced mode
    let grid = null;
    if (junctionConfig.mode === 'graded') {
        junctionConfig.nProfile.N0 = nConc;
        junctionConfig.pProfile.N0 = pConc;
        junctionConfig.nProfile.type = junctionConfig.profileType;
        junctionConfig.pProfile.type = junctionConfig.profileType;

        grid = solvePoisson({ Vbi, forwardBias, kT, mat, temperature });
    }

    return {
        // ... existing properties ...
        isGraded: junctionConfig.mode === 'graded',
        grid: grid,
        profileType: junctionConfig.profileType
    };
}
```

### Step 2.5: Modify updatePotentialChart()
**Location:** Line 614, add conditional logic

```javascript
function updatePotentialChart(calc) {
    let positions, phi_data, efield_Vcm, xp_um, xn_um, Vbi_eff;

    if (calc.isGraded && calc.grid) {
        // Use numerical solution
        const grid = calc.grid;
        positions = grid.x;
        phi_data = grid.phi;
        efield_Vcm = grid.E.map(e => e / 100);  // V/m → V/cm
        xp_um = grid.xp * 1e6;
        xn_um = grid.xn * 1e6;
        Vbi_eff = Math.max(...phi_data);

    } else {
        // Use analytical solution (existing code for abrupt junction)
        const dep = calculateDepletionParams(calc);
        // ... existing abrupt junction code (lines 615-663) ...
    }

    // Convert to micrometers
    const pos_um = positions.map(x => x * 1e6);

    // ... rest of chart rendering (existing code) ...
}
```

---

## PHASE 3: Band Diagram Integration (Lines 837-920)

### Step 3.1: Modify updateChart() to Use Numerical φ(x)
**Location:** Line 837

```javascript
function updateChart(calc) {
    const points = 100;
    const positions = [];
    const Ec_data = [], Ev_data = [], Ef_data = [], Ei_data = [];
    const Efn_data = [], Efp_data = [];

    if (calc.isGraded && calc.grid) {
        // Use numerical potential for band bending
        const grid = calc.grid;
        const N = grid.x.length;
        const step = Math.max(1, Math.floor(N / points));

        for (let i = 0; i < N; i += step) {
            const x_norm = grid.x[i] / (grid.xp + grid.xn) * 2;  // Normalize to [-1,1]
            positions.push(x_norm);

            // Band bending from electrostatic potential
            const qPhi = grid.phi[i];
            Ec_data.push(calc.Ec - qPhi);
            Ev_data.push(calc.Ev - qPhi);
            Ei_data.push(calc.Ei - qPhi);

            // Fermi level (constant in equilibrium)
            if (!calc.showQuasiFermi) {
                Ef_data.push(grid.x[i] < 0 ? calc.Ef_p : calc.Ef_n);
            } else {
                // Quasi-Fermi levels under bias
                const alpha = (grid.x[i] + grid.xp) / (grid.xp + grid.xn);
                Efn_data.push(calc.Efn_n * (1-alpha) + calc.Efn_p * alpha - qPhi);
                Efp_data.push(calc.Efp_n * (1-alpha) + calc.Efp_p * alpha - qPhi);
            }
        }

    } else {
        // Use existing ad-hoc bending for abrupt junction
        // ... existing code (lines 848-878) ...
    }

    // Chart rendering (existing code continues)
    // ... lines 879-1004 ...
}
```

---

## PHASE 4: Additional Profile Types (Lines 445-540)

### Step 4.1: Add Exponential Profile Support
Update `evaluateProfile()` function - already implemented in Phase 1 Step 1.2

### Step 4.2: Add Gaussian Profile Support
Update `evaluateProfile()` function - already implemented in Phase 1 Step 1.2

### Step 4.3: Add Custom Expression Support
```javascript
function evaluateCustomExpression(expr, x, N0) {
    // Safe evaluation of user-defined expressions
    try {
        const context = { x: x, N0: N0, Math: Math };
        const func = new Function(...Object.keys(context), `return ${expr};`);
        const result = func(...Object.values(context));

        // Validate
        if (isNaN(result) || !isFinite(result) || result < 0) {
            return N0;  // Fallback
        }
        return result;
    } catch (e) {
        console.error('Custom expression error:', e);
        return N0;
    }
}
```

Update `evaluateProfile()` to handle custom type:
```javascript
case 'custom':
    return evaluateCustomExpression(profile.expression, x_cm, profile.N0);
```

### Step 4.4: Add UI Controls for Each Profile Type
Add HTML elements and event listeners for exponential, gaussian, and custom parameters (similar structure to linear controls in Phase 2).

---

## PHASE 5: Visualization Enhancements (Lines 354-375, 1005-1280)

### Step 5.1: Add Doping Profile Chart
**Location:** After E-field chart (after line 353)

```html
<div id="dopingProfileSection" class="chart-container" style="display: none;">
    <h3>Doping Profile N_D(x), N_A(x)</h3>
    <div class="chart-wrapper" style="height: 400px;">
        <canvas id="dopingProfileChart"></canvas>
    </div>
</div>
```

### Step 5.2: Implement Doping Profile Chart Renderer
```javascript
function updateDopingProfileChart(calc) {
    if (!calc.isGraded || !calc.grid) {
        if (dopingProfileChart) dopingProfileChart.destroy();
        document.getElementById('dopingProfileSection').style.display = 'none';
        return;
    }

    document.getElementById('dopingProfileSection').style.display = 'block';

    const grid = calc.grid;
    const pos_um = grid.x.map(x => x * 1e6);
    const ND_cm3 = grid.ND.map(n => n / 1e6);
    const NA_cm3 = grid.NA.map(n => n / 1e6);

    if (dopingProfileChart) dopingProfileChart.destroy();

    const ctx = document.getElementById('dopingProfileChart').getContext('2d');
    dopingProfileChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: pos_um,
            datasets: [
                {
                    label: 'N_D(x)',
                    data: ND_cm3,
                    borderColor: '#2563eb',
                    backgroundColor: 'rgba(37,99,235,0.1)',
                    borderWidth: 2.5,
                    pointRadius: 0,
                    fill: true
                },
                {
                    label: 'N_A(x)',
                    data: NA_cm3,
                    borderColor: '#dc2626',
                    backgroundColor: 'rgba(220,38,38,0.1)',
                    borderWidth: 2.5,
                    pointRadius: 0,
                    fill: true
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                legend: { display: true, position: 'bottom' },
                annotation: {
                    annotations: createDepletionAnnotations({
                        xp: grid.xp, xn: grid.xn, W: grid.xp + grid.xn
                    })
                }
            },
            scales: {
                x: {
                    type: 'linear',
                    title: { display: true, text: 'Position (μm)' }
                },
                y: {
                    type: 'logarithmic',
                    title: { display: true, text: 'Doping (cm⁻³)' }
                }
            }
        }
    });
}
```

### Step 5.3: Update Main update() Function
**Location:** Line 1282

```javascript
function update() {
    const calc = calculate();
    updateInfo(calc);
    updateChart(calc);
    updatePotentialChart(calc);
    updateDopingProfileChart(calc);  // NEW
    updateCalculations(calc);

    document.getElementById('quasiFermiHelp').style.display =
        calc.showQuasiFermi ? 'list-item' : 'none';
}
```

---

## PHASE 6: Polish & Testing

### Step 6.1: Add Input Validation
```javascript
// Validate profile parameters
function validateProfileParams() {
    const type = junctionConfig.profileType;

    if (type === 'linear') {
        const aN = junctionConfig.nProfile.params.a;
        const aP = junctionConfig.pProfile.params.a;

        if (Math.abs(aN) < 1e15 || Math.abs(aN) > 1e25) {
            alert('Linear gradient out of range (1e15 - 1e25 cm⁻⁴)');
            return false;
        }
    }
    // Similar for exponential, gaussian

    return true;
}
```

### Step 6.2: Add Error Recovery
Already implemented in Phase 1 Step 1.6 (try-catch in `solvePoisson()`).

### Step 6.3: Testing Checklist

**Unit Tests (console verification):**
1. Linear profile with a=0 → matches abrupt junction
2. Exponential profile: verify doping decays correctly
3. Gaussian profile: verify peak at junction
4. Custom expression: test valid and invalid inputs
5. Thomas algorithm: verify solution satisfies A·φ = b

**Integration Tests (visual verification):**
1. Switch between abrupt/graded modes → no errors
2. Change profile type → charts update correctly
3. Adjust parameters → smooth updates
4. Edge cases: T=0K, T=700K, Va→Vbi

**Validation Against Theory:**
1. Linear graded: W ∝ V_bi^(1/3) (analytical formula exists)
2. Step function approximation → converges to abrupt results
3. Built-in potential: ∫E(x)dx = V_bi
4. Charge neutrality: ∫ρ(x)dx = 0

---

## Implementation Sequence (Critical Path)

### Minimum Viable Product (Phases 1-3)
1. ✅ Add data structures
2. ✅ Implement linear profile evaluation
3. ✅ Implement FDM Poisson solver
4. ✅ Add advanced mode checkbox
5. ✅ Add linear parameter sliders
6. ✅ Modify calculate() to call solver
7. ✅ Modify updatePotentialChart() for graded mode
8. ✅ Modify updateChart() for accurate band bending
9. ✅ Test linear graded junction end-to-end

**Estimated time:** 3-4 hours for MVP

### Full Feature Set (Phases 4-6)
10. Add exponential/gaussian/custom profiles (1 hour)
11. Add doping profile visualization (1 hour)
12. Add comprehensive error handling (30 min)
13. Testing and validation (1-2 hours)

**Total estimated time:** 6-8 hours for complete implementation

---

## Verification

After implementation, verify by:

1. **Open simulator in browser:** `open /Users/vblake/isd/pneg/pn-junction-simulator.html`

2. **Test abrupt mode (default):**
   - Verify all existing functionality works unchanged
   - Compare charts with previous version

3. **Enable Advanced Mode:**
   - Check "⚙️ Advanced: Graded Junctions"
   - Graded controls section appears

4. **Test linear profile:**
   - Set N_D = 1e16, N_A = 1e16 (base concentrations)
   - Adjust linear gradient sliders
   - Verify:
     - Doping profile chart shows linear variation
     - Potential φ(x) is smoother than abrupt case
     - E-field E(x) is non-linear (not triangular)
     - W increases compared to abrupt junction

5. **Test other profiles:**
   - Switch to Exponential → verify exponential decay
   - Switch to Gaussian → verify bell curve
   - Switch to Custom → test expression "N0 * (1 + 0.01*x)"

6. **Validate physics:**
   - Check ∫E(x)dx = V_bi (potential difference)
   - Verify E_c - E_v = E_g everywhere
   - Confirm E_F is constant in equilibrium

7. **Edge cases:**
   - Test at 0K, 180K, 300K, 700K
   - Test with forward bias approaching V_bi
   - Test with extreme doping asymmetry (N_D >> N_A)

8. **Browser compatibility:**
   - Test in Chrome, Firefox, Safari
   - Verify Charts.js annotation plugin loads
   - Check console for errors

---

## Risk Assessment

**Low Risk:**
- Backward compatibility preserved (abrupt mode unchanged)
- Fallback to abrupt mode on solver errors
- Pure additive feature (no breaking changes)

**Medium Risk:**
- Solver convergence for extreme parameters
- Performance with high grid resolution (500+ points)

**Mitigation:**
- Comprehensive error handling with user feedback
- Grid resolution can be adjusted in GRID_CONFIG
- Try-catch blocks prevent crashes

## Summary

This implementation adds full graded junction support to the P-N junction simulator while maintaining simplicity for beginners. The numerical Poisson solver enables accurate modeling of arbitrary doping profiles, and the advanced mode toggle keeps the interface clean.

**Core Innovation:**
- Replace analytical formulas with numerical FDM solver
- Support 4 profile types: linear, exponential, gaussian, custom
- Accurate band bending from electrostatic potential φ(x)

**User Benefits:**
- Explore realistic junction physics beyond textbook abrupt model
- Visualize how doping profiles affect E-field and depletion width
- Educational tool for semiconductor device physics courses
