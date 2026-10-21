"""N2O / ABS — RocketCEA が返せる値を全部出す

cea_obj_w_units.CEA_Obj の全 38 メソッド + fac_CR + 生出力を網羅する。
設計点は Pc=20 bar, O/F=6.5, eps=4。地上燃焼試験（海面上）を前提にサイジングする。
注意: hybrid_sim.py の現設計（Dt=12.1 mm, ox 0.200 kg/s）では燃焼中 Pc が
29.9 -> 28.8 bar になり、この Pc とは一致しない。
"""
import numpy as np
from rocketcea.cea_obj import add_new_fuel, add_new_oxidizer
from rocketcea.cea_obj_w_units import CEA_Obj

G0 = 9.80665

# ============================================================================
# 設計点  ← ここだけ書き換えれば良い
# ============================================================================
Pc, MR, eps = 20.0, 6.0, 4.0   # bar, -, -   
T_TANK    = 20.0               # °C   N2O タンク温度（自己加圧・飽和液）
P_AMB     = 1.01325            # bar  燃焼試験の外気圧（海面上）
ETA_CSTAR = 1.00               # -    c* 効率の仮定値。燃焼試験で実測して置き換える
F_TARGET  = 500.0              # N    海面上で欲しい推力

# ============================================================================
# 推進剤定義
# ============================================================================
add_new_fuel('ABS', """
fuel ABS  C 5.03  H 5.69  N 0.33  wt%=100.0
  h,cal=15507.0  t(k)=298.15  rho,g/cc=0.970
""")

# 飽和液 N2O の物性表（CoolProp 6.x の N2O 状態方程式から作成）
#   h は CEA 基準（RocketCEA の気体 'N2O' = 19467 cal/mol @298.15 K）に揃えた値。
#   同じ手順で標準沸点を計算すると 14528 cal/mol となり、
#   RocketCEA 組み込みの 'N2O_nbp'（14583 cal/mol）と 0.4 % 以内で一致する。
#   T[°C]   Psat[bar]  rho_l[kg/m3]  h[cal/mol]
N2O_SAT = np.array([
    [-20.0, 17.97, 998.1, 15831.8],
    [-10.0, 23.93, 954.9, 16044.8],
    [  0.0, 31.22, 907.1, 16269.6],
    [ 10.0, 40.01, 852.2, 16511.9],
    [ 15.0, 45.04, 820.7, 16642.7],
    [ 20.0, 50.53, 785.1, 16783.3],
    [ 25.0, 56.52, 742.9, 16939.0],
    [ 30.0, 63.08, 688.1, 17123.2],
])
Psat_tank, rho_ox, h_ox = (np.interp(T_TANK, N2O_SAT[:, 0], N2O_SAT[:, i]) for i in (1, 2, 3))
add_new_oxidizer('N2O_sat', f"""
oxid NitrousOxide  N 2.0 O 1.0  wt%=100.0
  h,cal={h_ox:.1f}  t(k)={T_TANK + 273.15:.2f}  rho,g/cc={rho_ox / 1000:.4f}
""")

SI = dict(isp_units='sec', cstar_units='m/s', pressure_units='Bar',
          temperature_units='K', sonic_velocity_units='m/s',
          enthalpy_units='kJ/kg', density_units='kg/m^3',
          specific_heat_units='kJ/kg-K', thermal_cond_units='W/cm-degC',
          viscosity_units='centipoise')

C = CEA_Obj(oxName='N2O_sat', fuelName='ABS', **SI)

hdr = lambda s: print(f'\n{"─"*70}\n {s}\n{"─"*70}')

print(f'╔{"═"*68}╗')
print(f'║ {C.get_description():<66} ║')
print(f'║ {f"Pc = {Pc:.0f} bar,  O/F = {MR:.1f},  eps = {eps:.0f},  N2O sat. liquid @ {T_TANK:.0f} C":<66} ║')
print(f'╚{"═"*68}╝')

# ============================================================================
hdr('0. 酸化剤の状態 — 標準沸点液 (N2O_nbp) との比較')
# ============================================================================
C_nbp = CEA_Obj(oxName='N2O_nbp', fuelName='ABS', **SI)
print(f'  タンク {T_TANK:.0f} °C の飽和液:  Psat = {Psat_tank:.2f} bar,  '
      f'rho = {rho_ox:.1f} kg/m3,  h = {h_ox:.0f} cal/mol')
print(f'  N2O_nbp (184.4 K)      :                      h = 14583 cal/mol')
print(f'  {"":24}{"Isp_vac[s]":>11}{"c*[m/s]":>10}{"Tc[K]":>9}')
for tag, obj in [(f'N2O 飽和液 {T_TANK:.0f} °C', C), ('N2O_nbp', C_nbp)]:
    i, c, t = obj.get_IvacCstrTc(Pc=Pc, MR=MR, eps=eps)
    print(f'  {tag:<24}{i:11.2f}{c:10.1f}{t:9.1f}')
print('  → 液が温かいぶんエンタルピーが高く、nbp を使うと性能を低めに見積もる')

# ============================================================================
hdr('1. 基本性能 — Isp / c* / Tc')
# ============================================================================
Isp_vac, cstar, Tc = C.get_IvacCstrTc(Pc=Pc, MR=MR, eps=eps)
print(f'  get_Isp            真空 Isp          {C.get_Isp(Pc=Pc, MR=MR, eps=eps):8.2f} s')
print(f'  get_Cstar          特性排気速度 c*   {C.get_Cstar(Pc=Pc, MR=MR):8.1f} m/s')
print(f'  get_Tcomb          燃焼温度          {C.get_Tcomb(Pc=Pc, MR=MR):8.1f} K')
print(f'  get_Throat_Isp     スロートまでのIsp {C.get_Throat_Isp(Pc=Pc, MR=MR):8.2f} s'
      f'   (eps=1 相当)')
print(f'  検算  Isp*g0 = c* x Cf …  {Isp_vac*G0:.1f} = {cstar:.1f} x '
      f'{Isp_vac*G0/cstar:.4f}')

# 平衡流 / 凍結流 / 凍結位置の3通り
Isp_fr_chm, _, _ = C.getFrozen_IvacCstrTc(Pc=Pc, MR=MR, eps=eps, frozenAtThroat=0)
Isp_fr_tht, _, _ = C.getFrozen_IvacCstrTc(Pc=Pc, MR=MR, eps=eps, frozenAtThroat=1)
print(f'\n  平衡流 (shifting)          {Isp_vac:8.2f} s   ← 理論上限')
print(f'  凍結流 (スロートで凍結)    {Isp_fr_tht:8.2f} s')
print(f'  凍結流 (燃焼室で凍結)      {Isp_fr_chm:8.2f} s   ← 理論下限')
d_frz = 1 - Isp_fr_chm / Isp_vac
print(f'  平衡-凍結の差 {100*d_frz:.1f} % は、c* 効率による損失'
      f' (1-eta = {100*(1-ETA_CSTAR):.0f} % 仮定) より{"小さい" if d_frz < 1-ETA_CSTAR else "大きい"}。')
print('  → 実機の見積もりは「平衡流 x c* 効率」で行い、凍結流の値を目安にしない')

# ============================================================================
hdr('2. 分子量と比熱比 — 3箇所')
# ============================================================================
for name, fn in [('chamber', C.get_Chamber_MolWt_gamma),
                 ('throat ', C.get_Throat_MolWt_gamma),
                 ('exit   ', C.get_exit_MolWt_gamma)]:
    mw, gam = fn(Pc=Pc, MR=MR, eps=eps)
    print(f'  {name}   MW = {mw:6.2f} g/mol,  gamma = {gam:.4f}')

# Isp と同時に取得する版（CEA 呼び出しが1回で済む）
print('\n  ワンショット版（同じ値を1回の CEA 実行で取る）:')
for label, fn in [('Chm', C.get_IvacCstrTc_ChmMwGam),
                  ('Tht', C.get_IvacCstrTc_ThtMwGam),
                  ('Exi', C.get_IvacCstrTc_exitMwGam)]:
    isp, cs, tc, mw, gam = fn(Pc=Pc, MR=MR, eps=eps)
    print(f'    {label}  Isp={isp:.2f}  c*={cs:.1f}  Tc={tc:.1f}  MW={mw:.2f}  γ={gam:.4f}')

# ============================================================================
hdr('3. 流れの状態量 — [chamber, throat, exit]')
# ============================================================================
rows = [
    ('温度        [K]     ', C.get_Temperatures(Pc=Pc, MR=MR, eps=eps), 1),
    ('密度        [kg/m3] ', C.get_Densities(Pc=Pc, MR=MR, eps=eps),    4),
    ('音速        [m/s]   ', C.get_SonicVelocities(Pc=Pc, MR=MR, eps=eps), 1),
    ('比熱 Cp     [kJ/kgK]', C.get_HeatCapacities(Pc=Pc, MR=MR, eps=eps), 4),
    ('エンタルピー[kJ/kg] ', C.get_Enthalpies(Pc=Pc, MR=MR, eps=eps),   1),
]
for label, vals, nd in rows:
    print(f'  {label}  ' + '  '.join(f'{v:>10.{nd}f}' for v in vals))
print('  ※ 比熱は平衡比熱（反応による寄与を含む）。伝熱計算には 6. の凍結比熱を使う')

print(f'\n  燃焼室だけを個別に取るショートカット:')
print(f'    get_Chamber_Density  {C.get_Chamber_Density(Pc=Pc, MR=MR, eps=eps):.4f} kg/m3')
print(f'    get_Chamber_Cp       {C.get_Chamber_Cp(Pc=Pc, MR=MR, eps=eps):.4f} kJ/kg-K (平衡)'
      f'   /  frozen=1: {C.get_Chamber_Cp(Pc=Pc, MR=MR, eps=eps, frozen=1):.4f}')
print(f'    get_Chamber_H        {C.get_Chamber_H(Pc=Pc, MR=MR, eps=eps):.1f} kJ/kg')
print(f'    get_Chamber_SonicVel {C.get_Chamber_SonicVel(Pc=Pc, MR=MR, eps=eps):.1f} m/s')

# ============================================================================
hdr('4. ノズル流れ — 圧力比・マッハ数')
# ============================================================================
PcOvPe = C.get_PcOvPe(Pc=Pc, MR=MR, eps=eps)
_, gam_c = C.get_Chamber_MolWt_gamma(Pc=Pc, MR=MR, eps=eps)
PcOvPt_ideal = ((gam_c + 1) / 2) ** (gam_c / (gam_c - 1))
print(f'  get_PcOvPe          Pc/Pe            {PcOvPe:8.2f}   -> Pe = {Pc/PcOvPe:.3f} bar')
print(f'  get_Throat_PcOvPe   Pc/Pt            {C.get_Throat_PcOvPe(Pc=Pc, MR=MR):8.4f}'
      f'   (完全気体・燃焼室 γ={gam_c:.3f} なら {PcOvPt_ideal:.4f})')
print(f'  get_MachNumber      出口マッハ数     {C.get_MachNumber(Pc=Pc, MR=MR, eps=eps):8.4f}')
print(f'  get_eps_at_PcOvPe   Pc/Pe=30 に必要な eps '
      f'{C.get_eps_at_PcOvPe(Pc=Pc, MR=MR, PcOvPe=30.0):6.3f}')

# ============================================================================
hdr('5. 外圧を考慮した性能 — 地上/真空/剥離判定')
# ============================================================================
# 注意: mode 文字列中の Pe / Psep は単位変換されず、常に psia で入っている
#       （RocketCEA 内部が psia 系のため。1 bar = 14.5038 psia）
for Pamb, tag in [(P_AMB, '海面上'), (0.5, '高度 ~5.5 km'), (1e-9, '真空')]:
    isp, mode = C.estimate_Ambient_Isp(Pc=Pc, MR=MR, eps=eps, Pamb=Pamb)
    cf_opt, cf_amb, _ = C.get_PambCf(Pamb=Pamb, Pc=Pc, MR=MR, eps=eps)
    print(f'  {tag:<12} Pamb={Pamb:9.5f} bar   Isp={isp:7.2f} s   Cf={cf_amb:.4f}   '
          f'[{mode}]  ← Pe は psia')

# get_PambCf / getFrozen_PambCf の戻り値は (CF_cea, CF_amb, mode)
#   CF_cea: CEA 自身の Cf = Pe=Pamb を仮定した値（圧力項を含まない）
#   CF_amb: 真空 Cf から外圧項 Pamb*eps/Pc を引いた、実際の外圧での Cf
#           （ただし Pe < Psep で剥離と判定されると Kalt-Badal の剥離モデル値になる）
cf_cea_eq, cf_amb_eq, _ = C.get_PambCf(Pamb=P_AMB, Pc=Pc, MR=MR, eps=eps)
cf_cea_fr, cf_amb_fr, _ = C.getFrozen_PambCf(Pamb=P_AMB, Pc=Pc, MR=MR, eps=eps)
print(f'\n  海面上の Cf           {"平衡":>8}{"凍結":>10}')
print(f'    外圧込み (CF_amb)   {cf_amb_eq:8.4f}{cf_amb_fr:10.4f}   ← 推力評価に使うのはこちら')
print(f'    CEA 値   (CF_cea)   {cf_cea_eq:8.4f}{cf_cea_fr:10.4f}   ← Pe=Pamb 仮定、圧力項なし')

# ============================================================================
hdr('6. 輸送物性 — 3箇所（アブレーティブ/ヒートシンク壁の熱負荷評価・Bartz 式に使う）')
# ============================================================================
# 小型ハイブリッドは再生冷却ではなくアブレーティブかヒートシンク構造が普通。
# Bartz 式などの対流熱伝達には凍結物性（frozen=1）を使うのが一般的。
print(f'  {"":8}  {"":6}{"Cp[kJ/kgK]":>11}  {"mu[cP]":>9}  {"k[W/cm-K]":>11}  {"Pr[-]":>7}')
for name, fn in [('chamber', C.get_Chamber_Transport),
                 ('throat ', C.get_Throat_Transport),
                 ('exit   ', C.get_Exit_Transport)]:
    for frz, ftag in [(1, '凍結'), (0, '平衡')]:
        cp, mu, k, pr = fn(Pc=Pc, MR=MR, eps=eps, frozen=frz)
        print(f'  {name}  {ftag:<4}  {cp:11.4f}  {mu:9.5f}  {k:11.6f}  {pr:7.4f}')

# ============================================================================
hdr('7. 化学組成 — 値は [injector, chamber, throat, exit] の4要素')
# ============================================================================
MIN_FRAC = 0.005
molWtD, moleD = C.get_SpeciesMoleFractions(Pc=Pc, MR=MR, eps=eps, min_fraction=MIN_FRAC)
_,      massD = C.get_SpeciesMassFractions(Pc=Pc, MR=MR, eps=eps, min_fraction=MIN_FRAC)
print(f'  {"species":<9}{"MW":>7}   {"mole@chm":>9}{"mole@exit":>10}   '
      f'{"mass@chm":>9}{"mass@exit":>10}')
# 注意: mole と mass では min_fraction の足切りが別々にかかるので
#       キーが一致しない（H2 のような軽い種は質量分率で落ちる）
fmt = lambda v: f'{v:.4f}' if v >= MIN_FRAC else '—'
for sp in sorted(moleD, key=lambda s: -moleD[s][1]):
    mf = massD.get(sp, [0.0] * 4)
    print(f'  {sp:<9}{molWtD[sp]:7.2f}   {fmt(moleD[sp][1]):>9}{fmt(moleD[sp][3]):>10}   '
          f'{fmt(mf[1]):>9}{fmt(mf[3]):>10}')
print(f'  — は足切り（< {MIN_FRAC}）で値が返らなかったもの。ゼロという意味ではない')

# ============================================================================
hdr('8. 当量比と混合比 — 最適 O/F の確認')
# ============================================================================
phi, r = C.get_eqratio(Pc=Pc, MR=MR, eps=eps)
print(f'  get_eqratio    O/F={MR} における当量比 phi = {phi:.4f}  (r = {r:.4f})')
print(f'                 phi > 1 は燃料リッチ')
for target in [1.0, 1.2, 1.5]:
    print(f'  getMRforER     phi = {target:.1f} にしたい -> O/F = '
          f'{C.getMRforER(ERphi=target):.4f}')
mr_stoich = C.getMRforER(ERphi=1.0)

# O/F スイープで最適値を確かめる（地上試験なので海面上 Isp を主指標にする）
mr_grid = np.arange(4.0, 9.01, 0.1)
isp_sl  = np.array([C.estimate_Ambient_Isp(Pc=Pc, MR=m, eps=eps, Pamb=P_AMB)[0] for m in mr_grid])
cs_grid = np.array([C.get_Cstar(Pc=Pc, MR=m) for m in mr_grid])
mr_opt_isp = mr_grid[isp_sl.argmax()]
mr_opt_cs  = mr_grid[cs_grid.argmax()]
print(f'\n  O/F スイープ ({mr_grid[0]:.1f}–{mr_grid[-1]:.1f}, 刻み 0.1):')
print(f'    海面上 Isp 最大  O/F = {mr_opt_isp:.1f}   Isp = {isp_sl.max():.2f} s')
print(f'    c* 最大          O/F = {mr_opt_cs:.1f}   c*  = {cs_grid.max():.1f} m/s')
isp_here = C.estimate_Ambient_Isp(Pc=Pc, MR=MR, eps=eps, Pamb=P_AMB)[0]
print(f'    設計点           O/F = {MR:.1f}   Isp = {isp_here:.2f} s  '
      f'(最大の {100*isp_here/isp_sl.max():.2f} %)')
print(f'\n  → 化学量論 O/F = {mr_stoich:.3f} に対し、設計 O/F = {MR} は'
      f' {MR/mr_stoich:.2f} 倍 = {"燃料リッチ" if phi > 1 else "酸化剤リッチ"}側'
      f'（phi = {phi:.2f}）')

# ============================================================================
hdr('9. 有限断面積燃焼室 (fac_CR) — 噴射器端と燃焼端の圧力差')
# ============================================================================
# fac_CR を指定すると、入力の Pc は「噴射器端」の圧力として扱われる（CEA の fac 仕様）。
# 燃焼端の圧力はレイリー流れ（加熱による全圧損失）で Pc より下がる。
CR = 6.0   # 収縮比 A_chamber / A_throat
Cfac = CEA_Obj(oxName='N2O_sat', fuelName='ABS', fac_CR=CR, **SI)
pinj_ratio = Cfac.get_Pinj_over_Pcomb(Pc=Pc, MR=MR)
print(f'  収縮比 CR = {CR}')
print(f'  get_Pinj_over_Pcomb   Pinj/Pcomb = {pinj_ratio:.5f}')
print(f'    -> 噴射器端 {Pc:.1f} bar なら燃焼端は {Pc/pinj_ratio:.2f} bar')
print(f'  get_Chamber_MachNumber  燃焼室マッハ数 = '
      f'{Cfac.get_Chamber_MachNumber(Pc=Pc, MR=MR):.5f}')
print(f'  Isp (fac_CR あり) = {Cfac.get_Isp(Pc=Pc, MR=MR, eps=eps):.2f} s'
      f'   / なし = {C.get_Isp(Pc=Pc, MR=MR, eps=eps):.2f} s')

# ============================================================================
hdr('10. 酸化剤タンク条件')
# ============================================================================
# ハイブリッドは燃料が燃焼室内のグレインなので、二液式のバルク密度・
# 密度比推力は使わない。タンクに関わるのは N2O の状態だけ。
print(f'  タンク温度 {T_TANK:.0f} °C  ->  蒸気圧 {Psat_tank:.2f} bar,  液密度 {rho_ox:.1f} kg/m3')
print(f'  噴射器差圧の余裕  Psat - Pc = {Psat_tank - Pc:.1f} bar'
      f'  (Pc の {100*(Psat_tank - Pc)/Pc:.0f} %)')
print('  ※ 自己加圧では放出中にタンクが冷えて Psat と密度が変わる')

# ============================================================================
hdr('11. 設計点のサイジング — 地上試験（海面上）')
# ============================================================================
#   F     = Cf_amb * Pc * At          （推力は外圧込みの Cf で評価する）
#   mdot  = Pc * At / (eta * c*_th)   （c* 効率を明示的に仮定）
#   → F と Pc を決めれば At は eta に依らず決まり、eta は必要流量に効く。
#     逆に酸化剤流量を固定すると、eta が下がるぶん Pc が下がる（hybrid_sim.py の形）。
Cf_sl  = cf_amb_eq
cstar_real = ETA_CSTAR * cstar
At     = F_TARGET / (Cf_sl * Pc * 1e5)
mdot   = Pc * 1e5 * At / cstar_real
mdot_ox, mdot_f = mdot * MR / (1 + MR), mdot / (1 + MR)
Isp_sl_real = F_TARGET / (mdot * G0)
Dt = 2 * np.sqrt(At / np.pi)
print(f'  海面上推力 {F_TARGET:.0f} N を Pc = {Pc:.0f} bar で出すなら:')
print(f'    Cf (海面上, 平衡)  {Cf_sl:.4f}   (真空 {Isp_vac*G0/cstar:.4f})')
print(f'    c* = {ETA_CSTAR:.2f} x {cstar:.1f} = {cstar_real:.1f} m/s  （c* 効率は仮定値）')
print(f'    Dt = {Dt*1000:.2f} mm,  De = {Dt*np.sqrt(eps)*1000:.2f} mm')
print(f'    総流量  {mdot:.4f} kg/s  (ox {mdot_ox:.4f} / fuel {mdot_f:.4f})')
print(f'    実効 Isp (海面上)  {Isp_sl_real:.1f} s')
print(f'    酸化剤体積流量  {mdot_ox/rho_ox*1000:.4f} L/s  (N2O {T_TANK:.0f} °C, {rho_ox:.1f} kg/m3)')
print('    fuel 流量はグレインの後退速度で決まる。上の値は O/F を保つための目標値')

print(f'\n  c* 効率に対する感度（F = {F_TARGET:.0f} N, Pc = {Pc:.0f} bar 固定）:')
for eta in [0.85, 0.90, 0.95, 1.00]:
    m_ = Pc * 1e5 * At / (eta * cstar)
    print(f'    eta = {eta:.2f}   mdot = {m_:.4f} kg/s   Isp_sl = {F_TARGET/(m_*G0):6.1f} s')

# ============================================================================
hdr('12. CEA の生出力（単位付きラッパには無いので .cea_obj 経由）')
# ============================================================================
txt = C.cea_obj.get_full_cea_output(Pc=Pc, MR=MR, eps=eps, subar=CR,
                                    pc_units='bar', output='siunits',
                                    short_output=0, show_transport=1,
                                    show_mass_frac=True)
open('n2o_abs_cea.txt', 'w').write(txt)
print(f'  {len(txt.splitlines())} 行を n2o_abs_cea.txt に保存')
print('  （subar=収縮比 を渡すと燃焼室断面の情報も出る）')

# ============================================================================
hdr('13. 開口比スイープ — eps と必要流量（F, Pc 固定・海面上）')
# ============================================================================
#   c* は eps に依らないので、mdot = F / (Cf_sl(eps) * eta * c*) となり
#   流量の変化は海面上 Cf の変化だけで決まる。Cf_sl は Pe = Pamb で最大 → 流量最小。
#   CF_amb は RocketCEA 内部で Pe < Psep（Kalt-Badal の剥離圧）になると剥離モデルに切り替わる。
#   その境界で曲線が折れる。剥離モデルは簡易なので、剥離域の値は目安程度。
import matplotlib.pyplot as plt

EPS_RANGE = (1.2, 10.0)   # スイープ範囲  ← ここを書き換える

eps_grid = np.linspace(*EPS_RANGE, 90)
_res     = [C.get_PambCf(Pamb=P_AMB, Pc=Pc, MR=MR, eps=e) for e in eps_grid]
cf_grid  = np.array([r[1] for r in _res])
sep_mask = np.array([r[2].startswith('Separated') for r in _res])
m_grid   = F_TARGET / (cf_grid * cstar_real)
eps_opt  = C.get_eps_at_PcOvPe(Pc=Pc, MR=MR, PcOvPe=Pc / P_AMB)
eps_sep  = eps_grid[sep_mask.argmax()] if sep_mask.any() else None   # 剥離モデルに入る最初の eps
m_at     = lambda e: np.interp(e, eps_grid, m_grid)

print(f'  適正膨張 (Pe = Pamb)        eps = {eps_opt:.3f}   mdot = {m_at(eps_opt):.4f} kg/s  (最小)')
if eps_sep is not None:
    pe_sep = Pc / C.get_PcOvPe(Pc=Pc, MR=MR, eps=eps_sep)
    print(f'  剥離モデル適用 (Kalt-Badal)  eps >= {eps_sep:.2f}  (Pe/Pamb = {pe_sep/P_AMB:.2f})')
else:
    print('  スイープ範囲内では剥離しない')
print(f'  設計点                      eps = {eps:.3f}   mdot = {mdot:.4f} kg/s'
      f'  (最小比 +{100*(mdot/m_at(eps_opt)-1):.2f} %)')
print(f'  {"eps":>6}{"Cf_sl":>9}{"mdot":>10}{"ox":>9}{"fuel":>9}{"Dt[mm]":>9}')
for e in [1.5, 2.0, 3.0, eps_opt, 4.0, 5.0, 6.0, 8.0, 10.0]:
    if not EPS_RANGE[0] <= e <= EPS_RANGE[1]:
        continue
    cf = C.get_PambCf(Pamb=P_AMB, Pc=Pc, MR=MR, eps=e)[1]
    m  = F_TARGET / (cf * cstar_real)
    dt = 2 * np.sqrt(F_TARGET / (cf * Pc * 1e5) / np.pi)
    print(f'  {e:6.2f}{cf:9.4f}{m:10.4f}{m*MR/(1+MR):9.4f}{m/(1+MR):9.4f}{dt*1000:9.2f}')

fig, ax = plt.subplots(figsize=(7.5, 5))
if eps_sep is not None:
    ax.axvspan(eps_sep, EPS_RANGE[1], color='#f3d9d2', alpha=0.6, lw=0)
    ax.text(eps_sep, 0.98, '  separated\n  (Kalt-Badal model)',
            va='top', fontsize=9, color='#a0442c', transform=ax.get_xaxis_transform())
ax.plot(eps_grid, m_grid,               color='#2a78d6', lw=2,   label='total')
ax.plot(eps_grid, m_grid * MR / (1 + MR), color='#eb6834', lw=1.6, label='oxidizer (N2O)')
ax.plot(eps_grid, m_grid / (1 + MR),      color='#1baf7a', lw=1.6, label='fuel (ABS)')
for e, lab, off in [(eps_opt, f'Pe = Pamb\nε = {eps_opt:.2f}', (0, 14)),
                    (eps,     f'design\nε = {eps:.2f}',        (0, -32))]:
    ax.plot(e, m_at(e), 'o', color='#2a78d6', ms=8, mec='white', mew=2, zorder=3)
    ax.annotate(f'{lab}\n{m_at(e):.4f} kg/s', xy=(e, m_at(e)), xytext=off,
                textcoords='offset points', ha='center', fontsize=9, color='#3d3d3a')
ax.set_xlabel('Expansion ratio  ε = Ae/At  [-]')
ax.set_ylabel('Mass flow rate  [kg/s]')
ax.set_title(f'Required mass flow vs ε   (F_sl = {F_TARGET:.0f} N, Pc = {Pc:.0f} bar, '
             f'O/F = {MR:.1f}, η_c* = {ETA_CSTAR:.2f})', fontsize=11, color='#1a1a19')
ax.set_xlim(*EPS_RANGE)
ax.set_ylim(bottom=0)
ax.grid(True, color='#e5e5e2', lw=0.8)
ax.set_axisbelow(True)
for side in ('top', 'right'):
    ax.spines[side].set_visible(False)
for side in ('left', 'bottom'):
    ax.spines[side].set_color('#c7c7c2')
ax.legend(loc='center right', frameon=False)
fig.tight_layout()
fig.savefig('mdot_vs_eps.png', dpi=150)
print('  グラフを mdot_vs_eps.png に保存')
plt.show()

# ============================================================================
hdr('14. 収縮比スイープ — CR と必要流量（F, 噴射器端 Pc, eps 固定・海面上）')
# ============================================================================
#   fac_CR ありの Pc は噴射器端圧力。CEA が返す c* は CR によらずほぼ一定だが、
#   これは実際の Pinj*At/mdot ではないので流量計算には使わない。
#   代わりにスロートの質量流束 G = rho_t * a_t（M=1）から直接求める:
#     F = mdot * Isp_vac * g0 - Pamb * eps * At,   mdot_ideal = G * At
#     -> At = F / (G * Isp_vac * g0 - Pamb * eps),   mdot = G * At / eta
#   CR -> ∞ で 11. のサイジングと一致する（剥離は考えない。eps=4 なら剥離域から遠い）。
CR_RANGE = (1.2, 10.0)   # スイープ範囲  ← ここを書き換える

def size_at(obj):
    """obj (CEA_Obj) について F_TARGET を出すための (mdot, At, Pcomb, Pe) を返す"""
    rho = obj.get_Densities(Pc=Pc, MR=MR, eps=eps)
    a   = obj.get_SonicVelocities(Pc=Pc, MR=MR, eps=eps)
    G   = rho[1] * a[1]                                   # kg/m2-s（スロート）
    isp = obj.get_Isp(Pc=Pc, MR=MR, eps=eps)
    At_ = F_TARGET / (G * isp * G0 - P_AMB * 1e5 * eps)
    pcomb = Pc / obj.get_Pinj_over_Pcomb(Pc=Pc, MR=MR) if obj.fac_CR else Pc
    pe  = Pc / obj.get_PcOvPe(Pc=Pc, MR=MR, eps=eps)
    return G * At_ / ETA_CSTAR, At_, pcomb, pe

cr_grid = np.linspace(*CR_RANGE, 70)
cr_res  = np.array([size_at(CEA_Obj(oxName='N2O_sat', fuelName='ABS', fac_CR=cr, **SI))
                    for cr in cr_grid])
m_inf, At_inf, _, _ = size_at(C)            # 無限断面積燃焼室（CR = ∞）
mcr, pcomb_cr = cr_res[:, 0], cr_res[:, 2]

print(f'  CR = ∞ (11. と同じ前提)   mdot = {m_inf:.4f} kg/s   Dt = {2*np.sqrt(At_inf/np.pi)*1000:.2f} mm')
print(f'  {"CR":>6}{"Pcomb[bar]":>12}{"Pcomb/Pinj":>12}{"Ma_chm":>9}{"mdot":>9}{"ox":>9}'
      f'{"fuel":>9}{"Dt[mm]":>9}{"Dc[mm]":>9}{"Δmdot":>9}')
for cr in [1.2, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0]:
    if not CR_RANGE[0] <= cr <= CR_RANGE[1]:
        continue
    o = CEA_Obj(oxName='N2O_sat', fuelName='ABS', fac_CR=cr, **SI)
    m, At_, pcb, pe = size_at(o)
    dt = 2 * np.sqrt(At_ / np.pi)
    print(f'  {cr:6.1f}{pcb:12.2f}{pcb/Pc:12.4f}{o.get_Chamber_MachNumber(Pc=Pc, MR=MR):9.4f}'
          f'{m:9.4f}{m*MR/(1+MR):9.4f}{m/(1+MR):9.4f}{dt*1000:9.2f}{dt*np.sqrt(cr)*1000:9.2f}'
          f'{100*(m/m_inf-1):+8.2f}%')
print('  ※ Dc はポート（燃焼室）内径の目安。ハイブリッドではポート径が燃焼中に広がるので、')
print('    実際の CR は燃焼開始時が最小で、時間とともに大きくなる。')

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True,
                               gridspec_kw=dict(height_ratios=[3, 2]))
ax1.plot(cr_grid, mcr, color='#2a78d6', lw=2)
ax1.axhline(m_inf, color='#2a78d6', lw=1, ls='--')
ax1.annotate(f'CR = ∞: {m_inf:.4f} kg/s', xy=(CR_RANGE[1], m_inf), xytext=(0, -4),
             textcoords='offset points', ha='right', va='top', fontsize=9, color='#3d3d3a')
for c_mark in [CR_RANGE[0], 2.0, CR]:
    m_ = np.interp(c_mark, cr_grid, mcr)
    ax1.plot(c_mark, m_, 'o', color='#2a78d6', ms=7, mec='white', mew=2, zorder=3)
    ax1.annotate(f'CR = {c_mark:.1f}\n{m_:.4f} kg/s (+{100*(m_/m_inf-1):.2f} %)',
                 xy=(c_mark, m_), xytext=(10, 2), textcoords='offset points',
                 va='bottom' if c_mark > CR_RANGE[0] else 'top', fontsize=9, color='#3d3d3a')
ax1.set_ylabel('Total mass flow rate  [kg/s]')
sec = ax1.secondary_yaxis('right', functions=(lambda m: 100 * (m / m_inf - 1),
                                             lambda d: m_inf * (1 + d / 100)))
sec.set_ylabel('vs CR = ∞  [%]', color='#6b6b66')
sec.spines['right'].set_color('#c7c7c2')
ax1.set_title(f'Required mass flow vs contraction ratio   (F_sl = {F_TARGET:.0f} N, '
              f'Pinj = {Pc:.0f} bar, O/F = {MR:.1f}, ε = {eps:.1f})', fontsize=10.5, color='#1a1a19')

ax2.plot(cr_grid, pcomb_cr, color='#8a5cc2', lw=2)
ax2.axhline(Pc, color='#8a5cc2', lw=1, ls='--')
ax2.annotate(f'Pinj = {Pc:.0f} bar', xy=(CR_RANGE[1], Pc), xytext=(0, -12),
             textcoords='offset points', ha='right', fontsize=9, color='#3d3d3a')
ax2.set_ylabel('Combustion-end pressure  [bar]')
ax2.set_xlabel('Contraction ratio  CR = Ac/At  [-]')
ax2.set_xlim(*CR_RANGE)
for ax in (ax1, ax2):
    ax.grid(True, color='#e5e5e2', lw=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color('#c7c7c2')
fig.tight_layout()
fig.savefig('mdot_vs_cr.png', dpi=150)
print('  グラフを mdot_vs_cr.png に保存')
plt.show()
