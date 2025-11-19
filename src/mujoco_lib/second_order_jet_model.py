import numpy as np


class SecondOrderJetModel:
    def __init__(self, dt, turbine_model=1):
        self._coeffs = np.array(
            [
                -1.78936773e-01,
                -4.02906083e00,
                -2.46787781e00,
                4.29405641e-01,
                -8.63414638e-01,
                -1.33204502e-02,
                3.80852208e00,
                9.40731803e-01,
                -1.48319385e00,
                -1.78761692e-01,
                1.62857455e-01,
                -1.23170475e00,
                4.64754367e-02,
            ]
        )
        self._mean_thrust = 107.907
        self._std_thrust = 71.388
        self._mean_throttle = 46.572
        self._std_throttle = 32.010
        self.set_turbine_model(turbine_model)
        self.dt = dt

    def set_turbine_model(self, turbine_model):
        if turbine_model == 1:
            self._mean_thrust = 107.907
            self._std_thrust = 71.388
        elif turbine_model == 2:
            self._mean_thrust = 123.35
            self._std_thrust = 83.58
        elif turbine_model == 3:
            self._mean_thrust = 143.90
            self._std_thrust = 97.51
        elif turbine_model == 4:
            self._mean_thrust = 226.13
            self._std_thrust = 153.23
        else:
            raise ValueError("Invalid turbine model. Choose from 1, 2, 3, or 4.")

    def set_turbine_coeffs(self, coeffs):
        if len(coeffs) != 13:
            raise ValueError("Coefficients array must have exactly 13 elements.")
        self._coeffs = np.array(coeffs)

    def thrust_noise_std(self, turbine_model):
        if turbine_model == 1:
            return 10.8
        elif turbine_model == 2:
            return 12.6
        elif turbine_model == 3:
            return 14.8
        elif turbine_model == 4:
            return 23.2
        else:
            raise ValueError("Invalid turbine model. Choose from 1, 2, 3, or 4.")

    def predict(self, T, T_dot, u):
        # T: current thrust (scalar)
        # T_dot: current thrust rate (scalar)
        # u: throttle command (scalar)

        # Standardize inputs
        T_std = self._standardize_thrust(T)
        T_dot_std = self._standardize_thrust_dot(T_dot)
        u_std = self._standardize_throttle(u)

        # Compute f, g, v
        f = self._compute_f(T_std, T_dot_std)
        g = self._compute_g(T_std, T_dot_std)
        v = self._compute_v(u_std)

        # Compute acceleration in standardized space
        T_dot_dot_std = f + g * v

        # Update standardized states
        T_dot_std_new = T_dot_std + T_dot_dot_std * self.dt
        T_std_new = T_std + T_dot_std_new * self.dt

        # Destandardize outputs
        T_new = self._destandardize_thrust(T_std_new)
        T_dot_new = self._destandardize_thrust_dot(T_dot_std_new)

        return T_new, T_dot_new

    def _standardize_thrust(self, T):
        return (T - self._mean_thrust) / self._std_thrust

    def _standardize_thrust_dot(self, T_dot):
        return T_dot / self._std_thrust

    def _standardize_throttle(self, u):
        return (u - self._mean_throttle) / self._std_throttle

    def _destandardize_thrust(self, T_std):
        return T_std * self._std_thrust + self._mean_thrust

    def _destandardize_thrust_dot(self, T_dot_std):
        return T_dot_std * self._std_thrust

    def _compute_f(self, thrust_std, thrust_dot_std):
        f = (
            self._coeffs[0]
            + self._coeffs[1] * thrust_std
            + self._coeffs[2] * thrust_dot_std
            + self._coeffs[3] * thrust_std * thrust_dot_std
            + self._coeffs[4] * thrust_std**2
            + self._coeffs[5] * thrust_dot_std**2
        )
        return f

    def _compute_g(self, thrust_std, thrust_dot_std):
        g = (
            self._coeffs[6]
            + self._coeffs[7] * thrust_std
            + self._coeffs[8] * thrust_dot_std
            + self._coeffs[9] * thrust_std * thrust_dot_std
            + self._coeffs[10] * thrust_std**2
            + self._coeffs[11] * thrust_dot_std**2
        )
        return g

    def _compute_v(self, throttle_std):
        v = throttle_std + self._coeffs[12] * throttle_std**2
        return v
