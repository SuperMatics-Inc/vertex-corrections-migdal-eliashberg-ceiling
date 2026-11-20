import Mathlib.Basic.Complex.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum

namespace VertexBound

-- In the basis (I, σ₁), conjugation by σ₃ negates the σ₁ coefficient.
-- The coefficient of t in the product of the three factors is the derivative.
theorem crossed_nambu (a b c t : ℂ) :
    (Complex.I * Complex.I * Complex.I + (-t*a)*(t*b)*Complex.I +
      (-t*a)*Complex.I*(-t*c) + Complex.I*(t*b)*(-t*c)) =
        -Complex.I + t^2 * (a*c-a*b-b*c) * Complex.I ∧
    ((-t*a)*Complex.I*Complex.I + Complex.I*(t*b)*Complex.I +
      Complex.I*Complex.I*(-t*c) + (-t*a)*(t*b)*(-t*c)) =
        t*(a-b+c) + t^3*a*b*c := by
  constructor
  · calc
      _ = Complex.I * Complex.I * Complex.I + t^2 * (a*c-a*b-b*c) * Complex.I := by ring
      _ = _ := by rw [Complex.I_mul_I, neg_one_mul]
  · calc
      _ = t*(-a+b-c)*(Complex.I * Complex.I) + t^3*a*b*c := by ring
      _ = _ := by rw [Complex.I_mul_I]; ring

theorem corner_average (a b c d : ℝ) :
    a*(a-b+c-d) + b*(b-a+d-c) + c*(c-d+a-b) + d*(d-c+b-a)
      = (a-b+c-d)^2 := by ring

theorem static_cancels (a b : ℝ) : a-a+b-b = 0 ∧ a-b+b-a = 0 := by
  constructor <;> ring

theorem four_square_schur (a b c d : ℝ) :
    (a-b+c-d)^2 ≤ 4*(a^2+b^2+c^2+d^2) := by
  nlinarith [sq_nonneg (a+b), sq_nonneg (a-c), sq_nonneg (a+d),
    sq_nonneg (b+c), sq_nonneg (b-d), sq_nonneg (c+d)]

theorem combine_errors (l q : ℝ) (hl : 0 ≤ l) (hq : 0 ≤ q) :
    q*(l/2+59/2) + (10/3)*q^2 ≤ 30*((1+l)*q+q^2) := by
  nlinarith [mul_nonneg hl hq, sq_nonneg q]

theorem scaled_guard (r e : ℝ) (hr : 0 < r) (he : 0 < e)
    (h : r*(1+e) ≤ 1) : r^2+e*r < 1 := by
  have hr' : r < 1 := by nlinarith [mul_pos hr he]
  nlinarith [mul_pos hr (sub_pos.mpr hr')]

-- Mutants: changing the middle insertion's sign leaves a spurious static term;
-- forgetting the corner average's 1/4 saturates the wrong norm at equal magnitudes.
theorem wrong_middle_survives : (1-1-1-1 : ℝ) ≠ 0 := by norm_num
theorem corner_factor_sharp : ((1 : ℝ)-(-1)+1-(-1))^2 = 16 := by norm_num

#print axioms crossed_nambu
#print axioms corner_average
#print axioms static_cancels
#print axioms four_square_schur
#print axioms combine_errors
#print axioms scaled_guard
#print axioms wrong_middle_survives
#print axioms corner_factor_sharp

end VertexBound
