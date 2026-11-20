import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.FieldSimp

namespace SignedFinite

-- G1-SIGNED/RIGOR.md:34,46: contraction and localized forward margins.
theorem forward_guards :
    (2*(1/16 : ℝ)+8*(1/16 : ℝ)^2 = 5/32) ∧
    (4*(1/16 : ℝ)+48*(1/16 : ℝ)^2 = 7/16) ∧
    (2*(2*(1/65536 : ℝ))+8*(2*(1/65536 : ℝ))^2 ≤ 5*(1/65536 : ℝ)) ∧
    ((1-5*(1/65536 : ℝ))/((1+5*(1/65536 : ℝ))^2+(8*(1/65536 : ℝ))^2) ≥ 1/2) ∧
    (8*(1/65536 : ℝ)/(1-5*(1/65536 : ℝ)) ≤ 16*(1/65536 : ℝ)) ∧
    (1/4096-(24576 : ℝ)/2^32-480/2^24 = 55/262144) ∧
    (55/262144 : ℝ) > 1/8192 := by norm_num

-- RIGOR.md:44: the witness boxes and the transfer window.
theorem forward_geometry (p q t : ℝ) (hp : -(1/16) ≤ p ∧ p ≤ 1/16)
    (hq : 1/2 ≤ q ∧ q ≤ 3/4) (ht : 1/2 ≤ t ∧ t ≤ 3/4) :
    p-q < -(1/16) ∧ p-t < -(1/16) ∧ p-q-t < -(1/16) ∧
    ((3/4 : ℝ)-1/2)/2 = 1/8 := by
  rcases hp with ⟨hp₁, hp₂⟩
  rcases hq with ⟨hq₁, hq₂⟩
  rcases ht with ⟨ht₁, ht₂⟩
  constructor
  · linarith
  constructor
  · linarith
  constructor
  · linarith
  · norm_num

-- RIGOR.md:42: all corners remain in the stated dispersion window.
theorem corner_window (p q t : ℝ) (hp : |p| ≤ 1/16)
    (hq : |q| ≤ 1) (ht : |t| ≤ 1) :
    |p-q| ≤ 2+1/16 ∧ |p-t| ≤ 2+1/16 ∧ |p-q-t| ≤ 2+1/16 := by
  rw [abs_le] at hp hq ht
  rcases hp with ⟨hp₁, hp₂⟩
  rcases hq with ⟨hq₁, hq₂⟩
  rcases ht with ⟨ht₁, ht₂⟩
  repeat' first | constructor | apply abs_le.mpr
  all_goals linarith

-- RIGOR.md:67: equal-argument integral coefficient and static guard.
theorem angular_guards :
    (3/8 : ℝ)-5/16 = 1/(2*(1+1)^3) ∧
    -(15/16 : ℝ)+250/1000 < -1/4 := by norm_num

-- RIGOR.md:40,67: witness temperature kernel and tail parameters.
theorem witness_guards :
    (1/(1+(2^13 : ℝ)^2) < 1/1000) ∧
    ((1/2^24 : ℝ) ≤ 1) ∧
    ((1827262477723295078215179/10^25 : ℝ) > 1/8) ∧
    ((1827262477723295078215179/10^25 : ℝ)*2^16 > 2^12) := by norm_num

-- RIGOR.md:67: absorption of the two-transfer angular error.
theorem angular_absorption (k s : ℝ) (hk : 0 < k) (hk' : k < 1/1000)
    (hs : 0 ≤ s) (hs' : s ≤ 5*k) :
    -(15/16)*k + 10*s^2 < -k/4 := by
  have hsq : s^2 ≤ (5*k)^2 := by nlinarith
  nlinarith

-- RIGOR.md:67: partial-fraction value used by the finite replay.
theorem angular_partial_fraction (a b : ℝ) (ha : a ≠ 0) (hb : b ≠ 0)
    (hd : b^2-a^2 ≠ 0) (hs : a+b ≠ 0) :
    let d := b^2-a^2
    let da := -a*b*a^2/d^2
    let db := -a*b*b^2/d^2
    let sa := (a*b-2*d*da)/d^2
    sa*(1/a-1/b)+da/(2*a^3)+db/(2*b^3) = 1/(2*(a+b)^3) := by
  dsimp
  field_simp
  ring

-- RIGOR.md:79: algebraic Jacobian, conditional on nonzero radii and width.
theorem overlap_jacobian (r rc d w : ℝ) (hr : r ≠ 0) (hc : rc ≠ 0)
    (hd : d ≠ 0) (hw : w ≠ 0) :
    2*w*3*r^2*((2*rc*d)/(2*r)) = 6*w*r*rc*d ∧
    (1/2 : ℝ)/(6*w*r*rc*d)/(2*w)^2 = 1/(48*w^3*r*rc*d) := by
  constructor <;> field_simp <;> ring

-- RIGOR.md:116: weak-coupling rational contraction endpoint.
theorem overlap_guards :
    (2*4096/(3*9) : ℝ) < 304 ∧
    (1/2^20 : ℝ)+4/2^10 < 1/128 ∧
    3*(1/128 : ℝ)+304*(1/128 : ℝ)^2 < 1/2 ∧
    4*(1/128 : ℝ)+6*304*(1/128 : ℝ)^2 < 1 := by norm_num

-- RIGOR.md:17,26: normal subtraction, before reanchoring.
def raw (rp ra rb rc qp qa qb qc yp ya yb yc : ℝ) : ℝ :=
  qp*yp*((qa*rb*rc+ra*qb*rc+ra*rb*qc-qa*qb*qc)*yp -
    ((rb*rc-qb*qc)*qa*ya+(ra*rc+qa*qc)*qb*yb+(ra*rb-qa*qb)*qc*yc))

-- RIGOR.md:22-26: exact algebraic four-anchor average, independent of integrals.
theorem four_anchor (rp ra rb rc qp qa qb qc yp ya yb yc : ℝ) :
    (raw rp ra rb rc qp qa qb qc yp ya yb yc +
     raw ra rp rc rb qa qp qc qb ya yp yc yb +
     raw rb rc rp ra qb qc qp qa yb yc yp ya +
     raw rc rb ra rp qc qb qa qp yc yb ya yp)/4 =
    -(qp*qa*qb*qc/4)*(yp-ya+yb-yc)^2 +
    (qp*qa*rb*rc*(yp-ya)^2+qp*qb*ra*rc*(yp-yb)^2+
     qp*qc*ra*rb*(yp-yc)^2 +
     qa*qb*rp*rc*(ya-yb)^2+qa*qc*rp*rb*(ya-yc)^2+
     qb*qc*rp*ra*(yb-yc)^2)/4 := by
  unfold raw
  ring

-- RIGOR.md:17: exact odd interpolation extracts a cubic polynomial's linear term.
theorem odd_interpolation (a b c d : ℝ) :
    (8*((a+b+c+d)-(a-b+c-d))-((a+2*b+4*c+8*d)-(a-2*b+4*c-8*d)))/12 = b := by ring

-- RIGOR.md:17: literal Nambu product, complex values as pairs.
def cmul (z w : ℝ × ℝ) : ℝ × ℝ :=
  (z.1*w.1-z.2*w.2, z.1*w.2+z.2*w.1)
def cadd (z w : ℝ × ℝ) : ℝ × ℝ := (z.1+w.1,z.2+w.2)
def mm (a b : Fin 2 → Fin 2 → ℝ × ℝ) : Fin 2 → Fin 2 → ℝ × ℝ :=
  fun i j => cadd (cmul (a i 0) (b 0 j)) (cmul (a i 1) (b 1 j))
def sigma : Fin 2 → Fin 2 → ℝ × ℝ :=
  fun i j => if i = 0 then (if j = 0 then (1,0) else (0,0))
             else (if j = 0 then (0,0) else (-1,0))
def line (r q y e : ℝ) : Fin 2 → Fin 2 → ℝ × ℝ :=
  fun i j => if i = 0 then (if j = 0 then (r,-q) else (-e*q*y,0))
             else (if j = 0 then (-e*q*y,0) else (-r,-q))
def vertex (ra rb rc qa qb qc ya yb yc e : ℝ) :=
  mm (mm (mm (mm (mm (mm sigma (line ra qa ya e)) sigma)
    (line rb qb yb e)) sigma) (line rc qc yc e)) sigma

theorem nambu_interpolation (ra rb rc qa qb qc ya yb yc : ℝ) :
    (8*((vertex ra rb rc qa qb qc ya yb yc 1 0 1).1 -
          (vertex ra rb rc qa qb qc ya yb yc (-1) 0 1).1) -
       ((vertex ra rb rc qa qb qc ya yb yc 2 0 1).1 -
          (vertex ra rb rc qa qb qc ya yb yc (-2) 0 1).1))/12 =
      (rb*rc-qb*qc)*qa*ya+(ra*rc+qa*qc)*qb*yb+
      (ra*rb-qa*qb)*qc*yc := by
  simp [vertex, mm, sigma, line, cadd, cmul, Fin.isValue]
  ring

-- RIGOR.md:26: the normal imaginary part of the same Nambu product.
theorem nambu_normal (ra rb rc qa qb qc ya yb yc : ℝ) :
    -(vertex ra rb rc qa qb qc ya yb yc 0 0 0).2 =
      qa*rb*rc+ra*qb*rc+ra*rb*qc-qa*qb*qc := by
  simp [vertex, mm, sigma, line, cadd, cmul, Fin.isValue]
  ring

#print axioms nambu_normal
#print axioms angular_absorption
#print axioms angular_partial_fraction
#print axioms nambu_interpolation
#print axioms forward_guards
#print axioms forward_geometry
#print axioms corner_window
#print axioms witness_guards
#print axioms angular_guards
#print axioms overlap_jacobian
#print axioms overlap_guards
#print axioms four_anchor
#print axioms odd_interpolation
end SignedFinite
