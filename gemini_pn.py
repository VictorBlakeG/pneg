

import sys
import numpy as np
import warnings

# Suppress warnings
warnings.filterwarnings("ignore")

from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QLabel, QLineEdit, QComboBox, QCheckBox, 
                             QGroupBox, QSizePolicy, QSlider)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
import matplotlib

# --- Material Constants (300K) ---
MATERIALS = {
    'Si': {'Eg': 1.12, 'eps': 11.7, 'ni': 1.0e10},
    'Ge': {'Eg': 0.66, 'eps': 16.0, 'ni': 2.4e13},
    'GaAs': {'Eg': 1.424, 'eps': 12.9, 'ni': 2.1e6},
    'In0.53Ga0.47As': {'Eg': 0.74, 'eps': 13.9, 'ni': 6.3e11},
    'InAs': {'Eg': 0.354, 'eps': 15.15, 'ni': 1.0e15}, 
    'GaSb': {'Eg': 0.726, 'eps': 15.7, 'ni': 1.5e12}
}

# Physical Constants
KB_EV = 8.617e-5   # eV/K
Q = 1.602e-19      # C
EPS0 = 8.854e-14   # F/cm

class PnJunctionGUI(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Interactive PN Junction: Final Version (TCAD Style)')
        self.showMaximized()

        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        self.main_layout = QVBoxLayout(self.main_widget)
        
        # --- Input Section ---
        self.inputs_group = QGroupBox("Simulation Parameters")
        self.inputs_group.setStyleSheet("font-size: 14pt; font-weight: bold;")
        self.inputs_layout = QHBoxLayout()
        self.inputs_group.setLayout(self.inputs_layout)
        
        # Col 1: Material
        self.col1_layout = QVBoxLayout()
        self.lbl_mat = QLabel("Material:")
        self.lbl_mat.setStyleSheet("font-weight: normal;")
        self.combo_mat = QComboBox()
        self.combo_mat.addItems(MATERIALS.keys())
        self.combo_mat.currentIndexChanged.connect(self.update_plots)
        self.col1_layout.addWidget(self.lbl_mat)
        self.col1_layout.addWidget(self.combo_mat)
        
        # Auto-Zoom Checkbox
        self.chk_zoom = QCheckBox("Auto-Zoom (3x W)")
        self.chk_zoom.setChecked(True)
        self.chk_zoom.setStyleSheet("font-weight: normal; font-size: 14pt; color: blue;")
        self.chk_zoom.toggled.connect(self.update_plots)
        self.col1_layout.addWidget(self.chk_zoom)
        
        self.inputs_layout.addLayout(self.col1_layout)

        # Col 2: Doping (Extended Range to 1e20)
        self.col2_layout = QVBoxLayout()
        # Range: 1e14 to 1e20 -> Log scale 14 to 20
        self.col2_layout.addWidget(self.create_slider_box("Na (cm^-3)", 14, 20, 15, self.update_na))
        self.col2_layout.addWidget(self.create_slider_box("Nd (cm^-3)", 14, 20, 15, self.update_nd))
        self.col2_layout.addWidget(self.create_slider_box("Diff Length (um)", 0.1, 10.0, 1.5, self.update_ldiff, step=0.1, is_log=False))
        self.inputs_layout.addLayout(self.col2_layout)

        # Col 3: Bias & Geometry
        self.col3_layout = QVBoxLayout()
        self.col3_layout.addWidget(self.create_slider_box("Applied Voltage (V)", -2.0, 1.0, 0.0, self.update_voltage, step=0.05, is_log=False))
        
        self.len_layout = QHBoxLayout()
        self.len_layout.addWidget(self.create_slider_box("P Length (um)", 0.5, 10.0, 3.0, self.update_lp, step=0.1, is_log=False))
        self.len_layout.addWidget(self.create_slider_box("N Length (um)", 0.5, 10.0, 3.0, self.update_ln, step=0.1, is_log=False))
        self.col3_layout.addLayout(self.len_layout)
        
        self.inputs_layout.addLayout(self.col3_layout)

        self.main_layout.addWidget(self.inputs_group)

        # --- Plotting Section ---
        matplotlib.rcParams.update({'font.size': 14}) # Updated to 14pt
        self.canvas = PlotCanvas(self, width=12, height=10)
        self.main_layout.addWidget(self.canvas)

        # Defaults
        self.Na = 1e15
        self.Nd = 1e15
        self.Va = 0.0
        self.Lp = 3.0
        self.Ln = 3.0
        self.L_diff = 1.5
        self.Temp = 300.0

        self.update_plots()

    def create_slider_box(self, label_text, min_val, max_val, default_val, callback, step=0.1, is_log=True):
        widget = QWidget()
        layout = QHBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        lbl = QLabel(label_text)
        lbl.setStyleSheet("font-weight: normal; font-size: 14pt;")
        layout.addWidget(lbl)
        txt = QLineEdit(str(default_val))
        txt.setFont(self.font()) # Inherit 14pt
        txt.setFixedWidth(110)
        layout.addWidget(txt)
        slider = QSlider(Qt.Horizontal)
        if is_log:
            slider.setRange(0, 100)
            val_norm = (np.log10(default_val) - min_val) / (max_val - min_val)
            slider.setValue(int(val_norm * 100))
        else:
            slider.setRange(int(min_val/step), int(max_val/step))
            slider.setValue(int(default_val/step))
        layout.addWidget(slider)

        def on_slider_change():
            if is_log:
                log_val = min_val + (slider.value()/100.0)*(max_val-min_val)
                real_val = 10**log_val
                if real_val > 1000: txt.setText(f"{real_val:.2e}")
                else: txt.setText(f"{real_val:.2f}")
            else:
                real_val = slider.value()*step
                txt.setText(f"{real_val:.2f}")
            callback(real_val)
        def on_text_change():
            try:
                val = float(txt.text())
                callback(val)
            except: pass
        slider.valueChanged.connect(on_slider_change)
        txt.returnPressed.connect(on_text_change)
        return widget

    def update_na(self, val): self.Na = val; self.update_plots()
    def update_nd(self, val): self.Nd = val; self.update_plots()
    def update_voltage(self, val): self.Va = val; self.update_plots()
    def update_lp(self, val): self.Lp = val; self.update_plots()
    def update_ln(self, val): self.Ln = val; self.update_plots()
    def update_ldiff(self, val): self.L_diff = val; self.update_plots()
    
    def update_plots(self):
        mat_name = self.combo_mat.currentText()
        props = MATERIALS[mat_name]
        is_zoom = self.chk_zoom.isChecked()
        self.canvas.plot_pn_junction(self.Na, self.Nd, self.Va, self.Lp, self.Ln, self.L_diff, self.Temp, props, is_zoom)


class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    def generate_smart_mesh(self, Lp_cm, Ln_cm, xp_cm, xn_cm, W_cm):
        # Dense Mesh for High Doping
        margin = 0.5 * W_cm # Wider margin for smoothness
        x_dense_min = -xp_cm - margin
        x_dense_max = xn_cm + margin
        
        x_dense_min = max(x_dense_min, -Lp_cm)
        x_dense_max = min(x_dense_max, Ln_cm)
        
        # Force 10,000 points in the active region
        N_dense = 10000
        x_dense = np.linspace(x_dense_min, x_dense_max, N_dense)
        
        x_bulk_p = np.array([])
        if x_dense_min > -Lp_cm:
            x_bulk_p = np.linspace(-Lp_cm, x_dense_min, 500, endpoint=False)
            
        x_bulk_n = np.array([])
        if x_dense_max < Ln_cm:
            x_bulk_n = np.linspace(x_dense_max, Ln_cm, 500)[1:] 
            
        x_combined = np.concatenate([x_bulk_p, x_dense, x_bulk_n])
        return np.sort(np.unique(x_combined))

    def solve_poisson_numerical(self, x_cm, Na, Nd, Va, ni, eps_s, T, Vbi):
        Vt = KB_EV * T
        N = len(x_cm)
        h = np.diff(x_cm)
        h_L = h[:-1]
        h_R = h[1:]
        avg_h = (h_L + h_R)
        
        N_dop = np.zeros(N)
        N_dop[x_cm <= 0] = -Na
        N_dop[x_cm > 0] = Nd

        phi_left = -Vt * np.log(Na/ni)
        phi_right = -Vt * np.log(Na/ni) + (Vbi - Va)

        phi = np.zeros(N)
        phi[x_cm <= 0] = phi_left
        phi[x_cm > 0] = phi_right
        
        center = np.searchsorted(x_cm, 0)
        sw = 200 # larger smooth window
        start = max(0, center-sw)
        end = min(N, center+sw)
        phi[start:end] = np.linspace(phi[start], phi[end-1], end-start)
        
        coeff = Q / eps_s
        phi_ref_h = phi_left + Vt*np.log(Na/ni) 
        phi_ref_e = phi_right - Vt*np.log(Nd/ni) 

        for _ in range(50): 
            p = ni * np.exp((phi_ref_h - phi)/Vt)
            n = ni * np.exp((phi - phi_ref_e)/Vt)
            rho = coeff * (p - n + N_dop)
            
            deriv_rho = coeff * (-p[1:-1]/Vt - n[1:-1]/Vt)
            
            a_vec = 2.0 / (h_L * avg_h)
            c_vec = 2.0 / (h_R * avg_h)
            b_fd  = -2.0 / (h_L * h_R)
            
            main_diag = b_fd + deriv_rho
            
            term_R = (phi[2:] - phi[1:-1]) / h_R
            term_L = (phi[1:-1] - phi[:-2]) / h_L
            d2phi = (2.0 / avg_h) * (term_R - term_L)
            F_interior = d2phi + rho[1:-1]
            
            M = N - 2
            c_prime = np.zeros(M)
            d_prime = np.zeros(M)
            
            c_prime[0] = c_vec[0] / main_diag[0]
            d_prime[0] = -F_interior[0] / main_diag[0]
            
            for i in range(1, M):
                temp = main_diag[i] - a_vec[i] * c_prime[i-1]
                if abs(temp) < 1e-30: temp = 1e-30
                c_prime[i] = c_vec[i] / temp
                d_prime[i] = (-F_interior[i] - a_vec[i] * d_prime[i-1]) / temp
                
            dphi = np.zeros(M)
            dphi[-1] = d_prime[-1]
            for i in range(M-2, -1, -1):
                dphi[i] = d_prime[i] - c_prime[i] * dphi[i+1]
            
            phi[1:-1] += dphi * 0.3 
            
            if np.max(np.abs(dphi)) < 1e-6: break
            
        return phi, Q*(p-n+N_dop), -np.gradient(phi, x_cm)

    def calculate_textbook_qfl(self, x_um, xp_um, xn_um, Va, Na, Nd, ni, T, E_i_curve, L_diff):
        Vt = KB_EV * T
        dE_p = -Vt * np.log(Na/ni)
        
        EFn = np.zeros_like(x_um)
        EFp = np.zeros_like(x_um)
        
        Ei_p_bulk = E_i_curve[0]
        EFp_bulk_val = Ei_p_bulk + dE_p 
        EFn_bulk_val = EFp_bulk_val + Va 
        
        for i, x in enumerate(x_um):
            if x < -xp_um:
                EFp[i] = EFp_bulk_val
                decay = np.exp((x + xp_um)/L_diff) 
                EFn[i] = EFp_bulk_val + (EFn_bulk_val - EFp_bulk_val) * decay
            elif x > xn_um:
                EFn[i] = EFn_bulk_val
                decay = np.exp(-(x - xn_um)/L_diff) 
                EFp[i] = EFn_bulk_val - (EFn_bulk_val - EFp_bulk_val) * decay
            else:
                EFn[i] = EFn_bulk_val
                EFp[i] = EFp_bulk_val
                
        return EFn, EFp

    def plot_pn_junction(self, Na, Nd, Va, Lp, Ln, L_diff, T, mat_props, is_zoom):
        self.fig.clear()
        
        gs = gridspec.GridSpec(3, 2, width_ratios=[1, 1.3], hspace=0.35, wspace=0.25)
        ax_rho = self.fig.add_subplot(gs[0, 0])
        ax_efield = self.fig.add_subplot(gs[1, 0], sharex=ax_rho)
        ax_pot = self.fig.add_subplot(gs[2, 0], sharex=ax_rho)
        ax_band = self.fig.add_subplot(gs[:, 1])

        ni = mat_props['ni']
        eps_s = mat_props['eps'] * EPS0
        Vt = KB_EV * T 
        Eg = mat_props['Eg']
        
        arg = (Na * Nd) / (ni**2)
        Vbi = Vt * np.log(arg) if arg > 1e-10 else 0.0
        
        Lp_cm = Lp * 1e-4
        Ln_cm = Ln * 1e-4
        
        eff_V = max(Vbi - Va, 1e-5) 
        W_approx = np.sqrt(2 * eps_s * eff_V * (1/Na + 1/Nd) / Q) 
        xp_approx = W_approx * Nd / (Na + Nd)
        xn_approx = W_approx * Na / (Na + Nd)
        
        x_cm = self.generate_smart_mesh(Lp_cm, Ln_cm, xp_approx, xn_approx, W_approx)
        x_um = x_cm * 1e4
        xp_um, xn_um = xp_approx*1e4, xn_approx*1e4

        phi_num, rho_real, E_real = self.solve_poisson_numerical(x_cm, Na, Nd, Va, ni, eps_s, T, Vbi)
        
        Ev_ref = 0.0
        V_shape = phi_num - phi_num[0]
        Ev_curve = Ev_ref - V_shape
        Ec_curve = Ev_curve + Eg
        Ei_curve = Ev_curve + Eg/2.0
        
        EFn_curve, EFp_curve = self.calculate_textbook_qfl(
            x_um, xp_um, xn_um, Va, Na, Nd, ni, T, Ei_curve, L_diff
        )

        rho_DA = np.zeros_like(x_um)
        rho_DA[(x_um > -xp_um) & (x_um <= 0)] = -Q * Na
        rho_DA[(x_um > 0) & (x_um < xn_um)] = Q * Nd
        
        ax_rho.plot(x_um, rho_DA, 'k--', alpha=0.5, label='DA')
        ax_rho.plot(x_um, rho_real, 'b-', lw=2, label='Num')
        ax_rho.set_ylabel(r'$\rho$ (C/cm$^3$)')
        ax_rho.tick_params(labelbottom=False)
        ax_rho.legend(fontsize=10)
        
        ax_efield.plot(x_um, E_real, 'r-', lw=2)
        ax_efield.set_ylabel('E-Field (V/cm)')
        ax_efield.tick_params(labelbottom=False)
        
        ax_pot.plot(x_um, V_shape, 'g-', lw=2)
        ax_pot.set_ylabel('Pot (V)')
        ax_pot.set_xlabel(r'x ($\mu m$)')

        ax_band.plot(x_um, Ec_curve, 'k-', lw=3, label=r'$E_c$')
        ax_band.plot(x_um, Ev_curve, 'k-', lw=3, label=r'$E_v$')
        ax_band.plot(x_um, Ei_curve, 'k--', lw=1.5, alpha=0.6, label=r'$E_i$')
        ax_band.plot(x_um, EFn_curve, 'c--', lw=2.5, label=r'$E_{Fn}$')
        ax_band.plot(x_um, EFp_curve, 'm--', lw=2.5, label=r'$E_{Fp}$')
        
        ax_band.axvspan(-xp_um, xn_um, color='gray', alpha=0.15)
        ax_band.axvline(-xp_um, color='gray', linestyle=':', linewidth=1)
        ax_band.axvline(xn_um, color='gray', linestyle=':', linewidth=1)
        
        y_vals = np.concatenate([Ec_curve, Ev_curve, EFn_curve, EFp_curve])
        y_min_plot, y_max_plot = np.min(y_vals), np.max(y_vals)
        y_span = y_max_plot - y_min_plot
        
        y_arrow = Ev_curve[0] + 0.15 * y_span 
        y_text = y_arrow + (0.05 * y_span)
        
        ax_band.annotate(
            "", xy=(-xp_um, y_arrow), xytext=(xn_um, y_arrow),
            arrowprops=dict(arrowstyle="<->", color='darkgreen', lw=2.0)
        )
        ax_band.text(0, y_text, f"$W={W_approx*1e4:.3f}\\mu m$",
                     color='darkgreen', ha='center',
                     fontsize=14, fontweight='bold')

        if abs(Va) > 0.05:
            mid = 0
            idx = np.abs(x_um - mid).argmin()
            y_n = EFn_curve[idx]
            y_p = EFp_curve[idx]
            ax_band.annotate("", xy=(0, y_n), xytext=(0, y_p),
                             arrowprops=dict(arrowstyle="<->", color='orange', lw=2))
            ax_band.text(0.1, (y_n+y_p)/2, f"qV={Va:.2f}eV", color='orange', fontweight='bold')

        ax_band.set_xlabel(r'Position ($\mu m$)')
        ax_band.set_ylabel('Energy (eV)')
        ax_band.set_title('Energy Band Diagram')
        ax_band.legend(loc='upper right', framealpha=0.9, fontsize=12)
        ax_band.grid(True, alpha=0.2)
        
        if is_zoom:
            view_w = max(W_approx*1e4 * 3.0, 0.2) # Scale 3x
            ax_band.set_xlim(-view_w, view_w)
            ax_rho.set_xlim(-view_w, view_w)
            ax_efield.set_xlim(-view_w, view_w)
            ax_pot.set_xlim(-view_w, view_w)
        else:
            ax_band.set_xlim(-Lp, Ln)
            ax_rho.set_xlim(-Lp, Ln)
            ax_efield.set_xlim(-Lp, Ln)
            ax_pot.set_xlim(-Lp, Ln)
            
        pad = 0.2
        ax_band.set_ylim(y_min_plot-pad, y_max_plot+pad)

        self.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyleSheet("QLabel{font-size: 14pt;} QLineEdit{font-size: 14pt;} QCheckBox{font-size: 14pt;}")
    gui = PnJunctionGUI()
    sys.exit(app.exec_())