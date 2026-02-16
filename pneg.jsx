jsx
import React, { useState, useMemo } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, Label } from 'recharts';

const MATERIALS = {
  Si: {
    name: 'Silicon',
    Eg0: 1.170,
    alpha: 4.73e-4,
    beta: 636,
    Nc: 2.8e19,
    Nv: 1.04e19,
    donors: ['P', 'As', 'Sb'],
    acceptors: ['B', 'Al', 'Ga', 'In'],
    donorEnergies: { P: 0.045, As: 0.054, Sb: 0.039 },
    acceptorEnergies: { B: 0.045, Al: 0.057, Ga: 0.065, In: 0.16 }
  },
  GaAs: {
    name: 'Gallium Arsenide',
    Eg0: 1.519,
    alpha: 5.41e-4,
    beta: 204,
    Nc: 4.7e17,
    Nv: 7.0e18,
    donors: ['Si', 'Ge', 'Sn', 'Se', 'Te'],
    acceptors: ['Si', 'Ge', 'Zn', 'Cd', 'Be', 'Mg'],
    donorEnergies: { Si: 0.0058, Ge: 0.006, Sn: 0.006, Se: 0.006, Te: 0.003 },
    acceptorEnergies: { Si: 0.035, Ge: 0.040, Zn: 0.031, Cd: 0.035, Be: 0.028, Mg: 0.028 }
  },
  Ge: {
    name: 'Germanium',
    Eg0: 0.7437,
    alpha: 4.77e-4,
    beta: 235,
    Nc: 1.04e19,
    Nv: 6.0e18,
    donors: ['P', 'As', 'Sb'],
    acceptors: ['B', 'Al', 'Ga', 'In'],
    donorEnergies: { P: 0.012, As: 0.013, Sb: 0.010 },
    acceptorEnergies: { B: 0.010, Al: 0.010, Ga: 0.011, In: 0.011 }
  },
  GaN: {
    name: 'Gallium Nitride',
    Eg0: 3.507,
    alpha: 9.14e-4,
    beta: 772,
    Nc: 2.3e18,
    Nv: 4.6e19,
    donors: ['Si', 'O', 'Ge'],
    acceptors: ['Mg', 'Zn', 'Cd', 'Be'],
    donorEnergies: { Si: 0.015, O: 0.029, Ge: 0.022 },
    acceptorEnergies: { Mg: 0.160, Zn: 0.340, Cd: 0.550, Be: 0.250 }
  },
  InP: {
    name: 'Indium Phosphide',
    Eg0: 1.421,
    alpha: 4.9e-4,
    beta: 327,
    Nc: 5.7e17,
    Nv: 1.1e19,
    donors: ['Si', 'S', 'Se', 'Sn', 'Te'],
    acceptors: ['Zn', 'Cd', 'Be', 'Mg'],
    donorEnergies: { Si: 0.007, S: 0.007, Se: 0.007, Sn: 0.007, Te: 0.007 },
    acceptorEnergies: { Zn: 0.040, Cd: 0.050, Be: 0.028, Mg: 0.040 }
  }
};

const PNJunctionSimulator = () => {
  const [material, setMaterial] = useState('Si');
  const [temperature, setTemperature] = useState(300);
  const [nDopant, setNDopant] = useState('P');
  const [pDopant, setPDopant] = useState('B');
  const [nConc, setNConc] = useState(1e16);
  const [pConc, setPConc] = useState(1e16);
  const [forwardBias, setForwardBias] = useState(0);
  const [showQuasiFermi, setShowQuasiFermi] = useState(false);

  const calculateBandgap = (T) => {
    const mat = MATERIALS[material];
    return mat.Eg0 - (mat.alpha * T * T) / (T + mat.beta);
  };

  const calculateNi = (T) => {
    const mat = MATERIALS[material];
    const Eg = calculateBandgap(T);
    const kT = 8.617e-5 * T;
    const Nc_T = mat.Nc * Math.pow(T / 300, 1.5);
    const Nv_T = mat.Nv * Math.pow(T / 300, 1.5);
    return Math.sqrt(Nc_T * Nv_T) * Math.exp(-Eg / (2 * kT));
  };

  const calculations = useMemo(() => {
    const mat = MATERIALS[material];
    const Eg = calculateBandgap(temperature);
    const kT = 8.617e-5 * temperature;
    const ni = calculateNi(temperature);
    
    const Ev = 0;
    const Ec = Eg;
    const Ei = Eg / 2;

    const Ef_n = Ec - kT * Math.log(mat.Nc * Math.pow(temperature / 300, 1.5) / nConc);
    const Ef_p = Ev + kT * Math.log(mat.Nv * Math.pow(temperature / 300, 1.5) / pConc);

    const Vbi = (Ef_n - Ef_p);

    const qV = forwardBias;
    const Efn_n = Ef_n;
    const Efp_p = Ef_p;
    const Efn_p = Ef_p + qV;
    const Efp_n = Ef_n - qV;

    return { Eg, Ev, Ec, Ei, Ef_n, Ef_p, Vbi, ni, Efn_n, Efp_p, Efn_p, Efp_n, kT };
  }, [material, temperature, nConc, pConc, forwardBias]);

  const bandData = useMemo(() => {
    const data = [];
    const points = 100;
    const { Ev, Ec, Ef_n, Ef_p, Vbi, Efn_n, Efp_p, Efn_p, Efp_n } = calculations;
    
    for (let i = 0; i <= points; i++) {
      const x = (i / points) * 2 - 1;
      let Ec_bent, Ev_bent, Ef, Efn, Efp;
      
      if (x < 0) {
        const bend = Math.abs(x) < 0.2 ? (1 - Math.abs(x) / 0.2) * Vbi / 2 : 0;
        Ec_bent = Ec - bend;
        Ev_bent = Ev - bend;
        Ef = Ef_n - bend;
        Efn = showQuasiFermi ? Efn_n - bend : null;
        Efp = showQuasiFermi ? Efp_n - bend : null;
      } else {
        const bend = Math.abs(x) < 0.2 ? (1 - Math.abs(x) / 0.2) * Vbi / 2 : Vbi / 2;
        Ec_bent = Ec - Vbi / 2 + bend;
        Ev_bent = Ev - Vbi / 2 + bend;
        Ef = Ef_p - Vbi / 2 + bend;
        Efn = showQuasiFermi ? Efn_p - Vbi / 2 + bend : null;
        Efp = showQuasiFermi ? Efp_p - Vbi / 2 + bend : null;
      }

      data.push({
        position: x,
        Ec: Ec_bent,
        Ev: Ev_bent,
        Ef: Ef,
        Efn: Efn,
        Efp: Efp,
        Ei: calculations.Ei - (x < 0 ? 0 : Vbi / 2) + (Math.abs(x) < 0.2 ? (1 - Math.abs(x) / 0.2) * Vbi / 2 : (x < 0 ? 0 : Vbi / 2))
      });
    }
    return data;
  }, [calculations, showQuasiFermi]);

  const mat = MATERIALS[material];

  return (

    <div className="w-full max-w-7xl mx-auto p-6 bg-gray-50">
      <h1 className="text-3xl font-bold mb-6 text-center text-gray-800">P-N Junction Fermi Level Simulator</h1>

      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 mb-6 bg-white p-6 rounded-lg shadow">
        {/* Column 1: Substrate, Temperature */}
        <div>
          <label className="block text-sm font-semibold mb-2 text-gray-700">Substrate Material</label>
          <select value={material} onChange={(e) => {
            setMaterial(e.target.value);
            setNDopant(MATERIALS[e.target.value].donors[0]);
            setPDopant(MATERIALS[e.target.value].acceptors[0]);
          }} className="w-full p-2 border rounded bg-white">
            {Object.keys(MATERIALS).map(m => (
              <option key={m} value={m}>{MATERIALS[m].name} ({m})</option>
            ))}
          </select>
          <div className="mt-4">
            <label className="block text-sm font-semibold mb-2 text-gray-700">Temperature: {temperature} K</label>
            <input type="range" min="0" max="700" value={temperature} onChange={(e) => setTemperature(parseInt(e.target.value))} className="w-full"/>
          </div>
        </div>

        {/* Column 2: P-type Dopant, ND */}
        <div>
          <label className="block text-sm font-semibold mb-2 text-gray-700">P-type Dopant</label>
          <select value={pDopant} onChange={(e) => setPDopant(e.target.value)} className="w-full p-2 border rounded bg-white">
            {mat.acceptors.map(a => (
              <option key={a} value={a}>{a} (Ea = {mat.acceptorEnergies[a]} eV)</option>
            ))}
          </select>
          <div className="mt-4">
            <label className="block text-sm font-semibold mb-2 text-gray-700">ND: {nConc.toExponential(2)} cm⁻³</label>
            <input type="range" min="14" max="19" step="0.1" value={Math.log10(nConc)} onChange={(e) => setNConc(Math.pow(10, parseFloat(e.target.value)))} className="w-full"/>
          </div>
        </div>

        {/* Column 3: N-type Dopant, NA */}
        <div>
          <label className="block text-sm font-semibold mb-2 text-gray-700">N-type Dopant</label>
          <select value={nDopant} onChange={(e) => setNDopant(e.target.value)} className="w-full p-2 border rounded bg-white">
            {mat.donors.map(d => (
              <option key={d} value={d}>{d} (Ed = {mat.donorEnergies[d]} eV)</option>
            ))}
          </select>
          <div className="mt-4">
            <label className="block text-sm font-semibold mb-2 text-gray-700">NA: {pConc.toExponential(2)} cm⁻³</label>
            <input type="range" min="14" max="19" step="0.1" value={Math.log10(pConc)} onChange={(e) => setPConc(Math.pow(10, parseFloat(e.target.value)))} className="w-full"/>
          </div>
        </div>

        {/* Forward Bias and Quasi-Fermi */}
        <div>
          <label className="block text-sm font-semibold mb-2 text-gray-700">Forward Bias: {forwardBias.toFixed(3)} V</label>
          <input type="range" min="0" max="3" step="0.01" value={forwardBias} onChange={(e) => setForwardBias(parseFloat(e.target.value))} className="w-full"/>
        </div>

        <div className="flex items-center">
          <label className="flex items-center cursor-pointer">
            <input type="checkbox" checked={showQuasiFermi} onChange={(e) => setShowQuasiFermi(e.target.checked)} className="mr-2 w-5 h-5"/>
            <span className="text-sm font-semibold text-gray-700">Show Quasi-Fermi Levels</span>
          </label>
        </div>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-4 mb-6">
        <div className="bg-white p-4 rounded-lg shadow">
          <h3 className="text-lg font-bold mb-3 text-gray-800">Varshni Parameters & Results</h3>
          <div className="space-y-2 text-sm">
            <div className="grid grid-cols-2 gap-2">
              <span className="font-semibold">Eg(0):</span><span>{mat.Eg0.toFixed(4)} eV</span>
              <span className="font-semibold">α:</span><span>{mat.alpha.toExponential(2)} eV/K</span>
              <span className="font-semibold">β:</span><span>{mat.beta} K</span>
              <span className="font-semibold text-blue-600">Eg({temperature}K):</span>
              <span className="text-blue-600 font-semibold">{calculations.Eg.toFixed(4)} eV</span>
              <span className="font-semibold">kT:</span><span>{calculations.kT.toFixed(5)} eV</span>
              <span className="font-semibold">ni:</span><span>{calculations.ni.toExponential(2)} cm⁻³</span>
            </div>
          </div>
        </div>

        <div className="bg-white p-4 rounded-lg shadow">
          <h3 className="text-lg font-bold mb-3 text-gray-800">Energy Levels (eV)</h3>
          <div className="space-y-2 text-sm">
            <div className="grid grid-cols-2 gap-2">
              <span className="font-semibold">Ec (Conduction):</span><span>{calculations.Ec.toFixed(4)} eV</span>
              <span className="font-semibold">Ev (Valence):</span><span>{calculations.Ev.toFixed(4)} eV</span>
              <span className="font-semibold">Ei (Intrinsic):</span><span>{calculations.Ei.toFixed(4)} eV</span>
              <span className="font-semibold text-green-600">EF,n (n-side):</span>
              <span className="text-green-600">{calculations.Ef_n.toFixed(4)} eV</span>
              <span className="font-semibold text-red-600">EF,p (p-side):</span>
              <span className="text-red-600">{calculations.Ef_p.toFixed(4)} eV</span>
              <span className="font-semibold">Vbi (Built-in):</span><span>{calculations.Vbi.toFixed(4)} V</span>
              {showQuasiFermi && (<>
                <span className="font-semibold text-purple-600">EFn split:</span>
                <span className="text-purple-600">{forwardBias.toFixed(4)} V</span>
              </>)}
            </div>
          </div>
        </div>
      </div>

      <div className="bg-white p-6 rounded-lg shadow">
        <h3 className="text-lg font-bold mb-4 text-gray-800 text-center">Energy Band Diagram</h3>
        <LineChart width={1000} height={500} data={bandData} margin={{ top: 20, right: 30, left: 60, bottom: 40 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="#ccc"/>
          <XAxis dataKey="position" type="number" domain={[-1, 1]} ticks={[-1, -0.5, 0, 0.5, 1]}
            tickFormatter={(val) => val === 0 ? 'Junction' : val < 0 ? 'n-side' : 'p-side'}>
            <Label value="Position" offset={-10} position="insideBottom"/>
          </XAxis>
          <YAxis domain={[Math.min(calculations.Ev, calculations.Ef_p - calculations.Vbi / 2) - 0.2, calculations.Ec + 0.2]}
            label={{ value: 'Energy (eV)', angle: -90, position: 'insideLeft' }}/>
          <Tooltip formatter={(value) => value?.toFixed(4) + ' eV'} labelFormatter={(label) => `Position: ${label.toFixed(2)}`}/>
          <Legend wrapperStyle={{ paddingTop: '20px' }}/>
          <Line type="monotone" dataKey="Ec" stroke="#2563eb" strokeWidth={3} name="Conduction Band (Ec)" dot={false}/>
          <Line type="monotone" dataKey="Ev" stroke="#dc2626" strokeWidth={3} name="Valence Band (Ev)" dot={false}/>
          <Line type="monotone" dataKey="Ei" stroke="#9333ea" strokeWidth={2} strokeDasharray="5 5" name="Intrinsic Level (Ei)" dot={false}/>
          {!showQuasiFermi && (<Line type="monotone" dataKey="Ef" stroke="#059669" strokeWidth={2} name="Fermi Level (Ef)" dot={false}/>)}
          {showQuasiFermi && (<>
            <Line type="monotone" dataKey="Efn" stroke="#0891b2" strokeWidth={2} strokeDasharray="3 3" name="Quasi-Fermi (electrons)" dot={false}/>
            <Line type="monotone" dataKey="Efp" stroke="#f59e0b" strokeWidth={2} strokeDasharray="3 3" name="Quasi-Fermi (holes)" dot={false}/>
          </>)}
        </LineChart>
        <div className="mt-4 text-sm text-gray-600 text-center">
          <p className="font-semibold">Material: {mat.name} | Temperature: {temperature} K | Bandgap: {calculations.Eg.toFixed(4)} eV</p>
          <p>N-type: {nDopant} ({nConc.toExponential(1)} cm⁻³) | P-type: {pDopant} ({pConc.toExponential(1)} cm⁻³)</p>
        </div>
      </div>

      <div className="mt-6 bg-blue-50 p-4 rounded-lg">
        <h4 className="font-bold text-gray-800 mb-2">Understanding the Diagram:</h4>
        <ul className="text-sm text-gray-700 space-y-1 list-disc list-inside">
          <li><strong>Blue line (Ec):</strong> Conduction band edge</li>
          <li><strong>Red line (Ev):</strong> Valence band edge</li>
          <li><strong>Purple dashed (Ei):</strong> Intrinsic Fermi level</li>
          <li><strong>Green line (EF):</strong> Equilibrium Fermi level</li>
          {showQuasiFermi && (<>
            <li><strong>Cyan dashed:</strong> Quasi-Fermi level for electrons (EFn)</li>
            <li><strong>Orange dashed:</strong> Quasi-Fermi level for holes (EFp)</li>
          </>)}
          <li><strong>Band bending:</strong> Shows the built-in electric field</li>
        </ul>
      </div>
    </div>
  );
};

export default PNJunctionSimulator;
