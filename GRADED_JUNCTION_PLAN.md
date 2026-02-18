# Graded Junction Implementation Plan

## Overview

Add spatially-varying doping profiles N_D(x) and N_A(x) to the P-N junction simulator to model realistic semiconductor junctions beyond the textbook abrupt model.

## Features

- **4 Profile Types**: Linear, Exponential, Gaussian, Custom expressions
- **Numerical Poisson Solver**: Full FDM solver for arbitrary doping profiles
- **Advanced Mode Toggle**: Preserves simple UI for beginners
- **Backward Compatible**: Existing abrupt junction mode unchanged

## Implementation Phases

### Phase 1: Core Numerical Solver (MVP)

**Location**: Lines 445-540

#### 1.1 Add Data Structures (after line 444)

```javascript
// Graded junction configuration
let junctionConfig = {
    mode: 'abrupt',  // 'abrupt' or 'graded'
    profileType: 'linear',
    nProfile: {
        type: 'constant',
        N0: null,
        params: { a: 1e20 }  // Linear gradient cm^-4
    },
    pProfile: {
        type: 'constant',
        N0: null,
        params: { a: -1e20 }
    }
};

const GRID_CONFIG = {
    points: 500,
    depletionFactor: 0.3
};

let gridData = {
    x: [], ND: [], NA: [], rho: [],
    phi: [], E: [], xp: 0, xn: 0
};
```

#### 1.2 Implement Profile Evaluation (after line 536)

```javascript
function evaluateProfile(profile, x_cm) {
    if (profile.type === 'constant') return profile.N0;

    const params = profile.params;
    switch (profile.type) {
        case 'linear':
            return Math.max(0, profile.N0 + params.a * x_cm);
        case 'exponential':
            return profile.N0 * Math.exp(params.lambda * x_cm);
        case 'gaussian':
            const arg = x_cm / params.sigma;
            return profile.N0 * Math.exp(-arg * arg);
        default:
            return profile.N0;
    }
}
```

#### 1.3 Implement Grid Creation

```javascript
function createAdaptiveGrid(ND_profile, NA_profile, Vbi, eps) {
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
        x, ND: new Array(N).fill(0), NA: new Array(N).fill(0),
        rho: new Array(N).fill(0), phi: new Array(N).fill(0),
        E: new Array(N).fill(0), xp: xp_est, xn: xn_est
    };
}
```

#### 1.4 Implement Poisson Solver (FDM)

```javascript
function solvePoissonFDM(x, rho, eps, Vbi_eff) {
    const N = x.length;
    const a = new Array(N).fill(0);
    const b = new Array(N).fill(0);
    const c = new Array(N).fill(0);
    const d = new Array(N).fill(0);

    // Boundary conditions
    b[0] = 1; c[0] = 0; d[0] = 0;
    b[N-1] = 1; a[N-1] = 0; d[N-1] = Vbi_eff;

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

#### 1.5 Implement E-field Calculation

```javascript
function computeElectricField(x, phi) {
    const N = x.length;
    const E = new Array(N);

    for (let i = 1; i < N - 1; i++) {
        E[i] = -(phi[i+1] - phi[i-1]) / (x[i+1] - x[i-1]);
    }

    E[0] = -(phi[1] - phi[0]) / (x[1] - x[0]);
    E[N-1] = -(phi[N-1] - phi[N-2]) / (x[N-1] - x[N-2]);

    return E;
}
```

#### 1.6 Implement Main Solver

```javascript
function solvePoisson(calc) {
    const ND_profile = junctionConfig.nProfile;
    const NA_profile = junctionConfig.pProfile;
    const eps = eps0 * eps_r[currentMaterial];
    const Vbi_eff = Math.max(calc.Vbi - calc.forwardBias, 0.001);

    try {
        const grid = createAdaptiveGrid(ND_profile, NA_profile, Vbi_eff, eps);

        for (let i = 0; i < grid.x.length; i++) {
            const x_cm = grid.x[i] * 100;
            grid.ND[i] = evaluateProfile(ND_profile, x_cm) * 1e6;
            grid.NA[i] = evaluateProfile(NA_profile, x_cm) * 1e6;
            grid.rho[i] = q * (grid.ND[i] - grid.NA[i]);
        }

        grid.phi = solvePoissonFDM(grid.x, grid.rho, eps, Vbi_eff);
        grid.E = computeElectricField(grid.x, grid.phi);

        // Find depletion edges
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

### Phase 2: UI Integration

#### 2.1 Add Advanced Mode Checkbox (after line 302)

```html
<div class="control-group checkbox-group">
    <input type="checkbox" id="advancedMode">
    <label for="advancedMode">⚙️ Advanced: Graded Junctions</label>
</div>
```

#### 2.2 Add Graded Controls Section (after line 321)

```html
<div id="gradedControls" class="controls" style="display: none; background: #fef3c7; border: 2px solid #f59e0b;">
    <h4 style="grid-column: 1 / -1; color: #92400e;">📊 Graded Junction Configuration</h4>

    <div class="control-group">
        <label for="profileType">Profile Type</label>
        <select id="profileType">
            <option value="linear">Linear: N(x) = N₀ + a·x</option>
            <option value="exponential">Exponential: N(x) = N₀·e^(λx)</option>
            <option value="gaussian">Gaussian: N(x) = N₀·e^(-(x/σ)²)</option>
            <option value="custom">Custom Expression</option>
        </select>
    </div>

    <div id="linearParamsN" class="control-group">
        <label>N-side gradient (a): <span id="linearNValue">1.00e+20</span> cm⁻⁴</label>
        <input type="range" id="linearN" min="18" max="22" step="0.1" value="20">
    </div>

    <div id="linearParamsP" class="control-group">
        <label>P-side gradient (a): <span id="linearPValue">-1.00e+20</span> cm⁻⁴</label>
        <input type="range" id="linearP" min="-22" max="-18" step="0.1" value="-20">
    </div>

    <div class="control-group" style="grid-column: 1 / -1;">
        <div style="background: #dbeafe; padding: 10px; border-radius: 4px; font-size: 0.875rem;">
            <strong>Note:</strong> N_D and N_A sliders set base concentration N₀.
        </div>
    </div>
</div>
```

#### 2.3 Add Event Listeners (after line 1330)

```javascript
document.getElementById('advancedMode').addEventListener('change', function() {
    junctionConfig.mode = this.checked ? 'graded' : 'abrupt';
    document.getElementById('gradedControls').style.display = this.checked ? 'grid' : 'none';
    update();
});

document.getElementById('profileType').addEventListener('change', function() {
    junctionConfig.profileType = this.value;
    ['linearParamsN', 'linearParamsP'].forEach(id => {
        document.getElementById(id).style.display = this.value === 'linear' ? 'flex' : 'none';
    });
    update();
});

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

#### 2.4 Modify calculate() (line 483)

Add before return statement:

```javascript
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
```

#### 2.5 Modify updatePotentialChart() (line 614)

Add at beginning of function:

```javascript
let positions, phi_data, efield_Vcm, xp_um, xn_um, Vbi_eff;

if (calc.isGraded && calc.grid) {
    const grid = calc.grid;
    positions = grid.x;
    phi_data = grid.phi;
    efield_Vcm = grid.E.map(e => e / 100);
    xp_um = grid.xp * 1e6;
    xn_um = grid.xn * 1e6;
    Vbi_eff = Math.max(...phi_data);
} else {
    const dep = calculateDepletionParams(calc);
    // ... existing code for abrupt junction ...
}
```

---

### Phase 3: Band Diagram Integration

#### 3.1 Modify updateChart() (line 837)

Replace band bending calculation:

```javascript
if (calc.isGraded && calc.grid) {
    const grid = calc.grid;
    const N = grid.x.length;
    const step = Math.max(1, Math.floor(N / points));

    for (let i = 0; i < N; i += step) {
        const x_norm = grid.x[i] / (grid.xp + grid.xn) * 2;
        positions.push(x_norm);

        const qPhi = grid.phi[i];
        Ec_data.push(calc.Ec - qPhi);
        Ev_data.push(calc.Ev - qPhi);
        Ei_data.push(calc.Ei - qPhi);

        if (!calc.showQuasiFermi) {
            Ef_data.push(grid.x[i] < 0 ? calc.Ef_p : calc.Ef_n);
        } else {
            const alpha = (grid.x[i] + grid.xp) / (grid.xp + grid.xn);
            Efn_data.push(calc.Efn_n * (1-alpha) + calc.Efn_p * alpha - qPhi);
            Efp_data.push(calc.Efp_n * (1-alpha) + calc.Efp_p * alpha - qPhi);
        }
    }
} else {
    // ... existing ad-hoc bending code ...
}
```

---

### Phase 4: Additional Profile Types

#### 4.1 Add Custom Expression Support

```javascript
function evaluateCustomExpression(expr, x, N0) {
    try {
        const context = { x: x, N0: N0, Math: Math };
        const func = new Function(...Object.keys(context), `return ${expr};`);
        const result = func(...Object.values(context));

        if (isNaN(result) || !isFinite(result) || result < 0) {
            return N0;
        }
        return result;
    } catch (e) {
        console.error('Custom expression error:', e);
        return N0;
    }
}
```

Update evaluateProfile() to add:

```javascript
case 'custom':
    return evaluateCustomExpression(profile.expression, x_cm, profile.N0);
```

#### 4.2 Add UI for Exponential/Gaussian/Custom

Add similar control groups for each profile type with appropriate parameter ranges.

---

### Phase 5: Visualization Enhancements

#### 5.1 Add Doping Profile Chart (after line 353)

```html
<div id="dopingProfileSection" class="chart-container" style="display: none;">
    <h3>Doping Profile N_D(x), N_A(x)</h3>
    <div class="chart-wrapper" style="height: 400px;">
        <canvas id="dopingProfileChart"></canvas>
    </div>
</div>
```

#### 5.2 Implement Chart Renderer

```javascript
let dopingProfileChart = null;

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

#### 5.3 Update main update() (line 1282)

Add:

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

### Phase 6: Polish & Testing

#### 6.1 Testing Checklist

**Unit Tests:**
- Linear profile with a=0 matches abrupt junction
- Exponential profile decays correctly
- Gaussian profile peaks at junction
- Custom expressions validate properly

**Integration Tests:**
- Switch between modes without errors
- All profile types update charts correctly
- Temperature and bias changes work in graded mode

**Physics Validation:**
- ∫E(x)dx = V_bi
- E_c - E_v = E_g everywhere
- E_F constant in equilibrium

---

## Quick Start

1. Open `/Users/vblake/isd/pneg/pn-junction-simulator.html`
2. Implement Phase 1 (solver foundation)
3. Test solver in console
4. Implement Phase 2 (UI)
5. Test linear profile visually
6. Continue through remaining phases

## Estimated Time

- **MVP (Phases 1-3):** 3-4 hours
- **Full (Phases 4-6):** 6-8 hours total

## Key Files

- **Main file:** `/Users/vblake/isd/pneg/pn-junction-simulator.html`
- **Plan:** `/Users/vblake/.claude/plans/starry-dreaming-scroll.md`

## Success Criteria

✅ Advanced mode checkbox toggles graded controls
✅ Linear profile shows non-triangular E-field
✅ All 4 profile types work
✅ Doping chart displays
✅ Abrupt mode still works perfectly
✅ No console errors
