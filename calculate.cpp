#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <sstream>
#include <functional>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Convert degrees to radians
static double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}
// Convert radians to degrees
static double rad2deg(double rad) {
    return rad * 180.0 / M_PI;
}

// Ensure theta within [-90, 90]
static double constrain_theta(double theta) {
    if (theta > 90.0) return 90.0;
    if (theta < -90.0) return -90.0;
    return theta;
}

// Effective k as function of theta (deg)
static double calculate_k_effective(double theta_deg) {
    const double k_base = 0.002138;
    const double c1 = 8.308074e-4;
    const double c2 = 9.296451e-5;
    return k_base * (1 + c1 * theta_deg + c2 * theta_deg * theta_deg);
}

// Solve f(theta)=0 using the secant method
static bool solve_theta(const std::function<double(double)>& func,
    double x0, double x1, double& root) {
    const int max_iter = 50;
    const double tol = 1e-6;
    double f0 = func(x0);
    double f1 = func(x1);
    for (int i = 0; i < max_iter; ++i) {
        if (std::fabs(f1 - f0) < 1e-12) break;
        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        double f2 = func(x2);
        if (std::fabs(f2) < tol) {
            root = x2;
            return true;
        }
        x0 = x1; f0 = f1;
        x1 = x2; f1 = f2;
    }
    root = x1;
    return false;
}

// Calculate original theta given r, alpha_rad, alpha_deg and initial guess (deg)
static double calculate_original_theta(double r, double alpha_rad,
    double alpha_deg, double init_theta_deg) {
    double k_eff = calculate_k_effective(init_theta_deg);
    auto func = [&](double theta_rad) {
        double s = std::sin(theta_rad);
        return 2 * std::sin(alpha_rad) * (1 - s * s)
            - std::cos(alpha_rad) * (std::sin(2 * theta_rad) - k_eff * r * std::cos(alpha_rad));
        };
    // Initial guesses around init_theta_deg
    double t0 = deg2rad(init_theta_deg - 5.0);
    double t1 = deg2rad(init_theta_deg + 5.0);
    double theta_rad;
    if (!solve_theta(func, t0, t1, theta_rad)) {
        // fallback to guess
        return init_theta_deg;
    }
    return rad2deg(theta_rad);
}

// Low-theta correction
static double calculate_correction_low(const std::vector<double>& p,
    double r, double a) {
    double corr = 0.0;
    corr += p[0] * r + p[1] * a * a + p[2] * std::pow(a, 3) + p[3] * std::pow(a, 4);
    corr += p[4] * r * r + p[5] * std::pow(r, 3) + p[6] * r * a + p[7] * r * r * a;
    corr += p[8] * r * a * a + p[9] * std::sin(a) + p[10] * std::cos(a) + p[11] * std::fabs(a);
    corr += p[12] * std::exp(-r / 1000.0) + p[13] * std::exp(-r / 500.0) + p[14] * std::exp(-r / 250.0);
    corr += p[15] * std::exp(-std::fabs(a) / 10.0) + p[16] * std::exp(-std::fabs(a) / 5.0);
    corr += p[17] * std::log(r + 1.0) + p[18] * std::log(std::fabs(a) + 1.0);
    corr += p[19] * std::sin(r / 100.0) + p[20] * std::cos(r / 100.0);
    corr += p[21] * std::sin(a / 10.0) + p[22] * std::cos(a / 10.0);
    corr += p[23] * std::tanh(r / 300.0) + p[24] * std::tanh(a / 5.0);
    corr += p[25] * r * std::sin(a / 10.0) + p[26] * r * std::cos(a / 10.0);
    corr += p[27] * std::sqrt(r) + p[28] * std::sqrt(std::fabs(a));
    corr += p[29] * std::exp(-r / 150.0) + p[30] * std::exp(-r / 100.0);
    if (r < 150) {
        corr += p[31] * (150 - r) / 150.0
            + p[32] * std::sin(M_PI * r / 150.0)
            + p[33] * std::cos(M_PI * r / 150.0)
            + p[34] * a * a
            + p[35] * r * a
            + p[36] * std::exp(-std::fabs(a) / 3.0)
            + p[37] * std::log(r / 150.0 + 0.1);
    }
    // clamp to [-30, 30]
    corr = std::max(-30.0, std::min(30.0, corr));
    return corr;
}

// High-theta correction
static double calculate_correction_high(const std::vector<double>& p,
    double r, double a) {
    double corr = 0.0;
    corr += p[0] * r + p[1] * a * a + p[2] * std::pow(a, 3) + p[3] * r * r + p[4] * r * a;
    corr += p[5] * std::sin(deg2rad(a)) + p[6] * std::cos(deg2rad(a));
    corr += p[7] * std::sin(deg2rad(2 * a)) + p[8] * std::cos(deg2rad(2 * a));
    corr += p[9] * std::exp(-r / 100.0) + p[10] * std::exp(-r / 200.0) + p[11] * std::exp(-r / 300.0);
    corr += p[12] * std::exp(-std::fabs(a) / 5.0) + p[13] * std::exp(-std::fabs(a) / 10.0) + p[14] * std::exp(-std::fabs(a) / 15.0);
    corr += p[15] * r * std::sin(deg2rad(a / 10.0)) + p[16] * r * std::cos(deg2rad(a / 10.0));
    corr += p[17] * std::sqrt(r) * std::tanh(a / 10.0);
    corr += p[18] * ((a > 45.0) ? std::pow((a - 45.0) / 45.0, 2) : 0.0)
        + p[19] * ((r > 200.0) ? (r - 200.0) / 100.0 : 0.0);
    // clamp to [-20, 20]
    corr = std::max(-20.0, std::min(20.0, corr));
    return corr;
}

int main() {
    std::vector<double> low_params = {
        0.05906051, 0.01922017, -0.09128481, 0.00227226, -0.00034493,
        0.00000029, 0.00505409, -0.00009850, 0.00455790, -0.05668212,
        0.01808962, 0.02486718, 0.00000941, 0.01426906, 0.01337904,
        -0.00347473, 0.00585867, -0.00318687, 0.00620625, -0.01442806,
        0.00265681, 0.01354817, 0.01208332, 0.01058386, -0.00150257,
        0.01667214, 0.00983190, -0.00946299, 0.03288948, 0.00304788,
        0.00352337, 0.00714029, -0.00714395, -0.03475796, 0.01317679,
        0.08707241, 0.03757451, 0.00491534
    };
    std::vector<double> high_params = {
        0.31953025, -0.00230637, 0.00030090, 0.00005720, -0.00210798,
        -0.15638930, 0.72346818, -0.19508769, 0.69479213, 0.43230354,
        0.67890975, 0.60368317, 0.83926671, -1.21558890, -1.10351086,
        0.22915913, -0.33235095, 0.26818104, 0.08942329, 1.45841379
    };

    std::string line;
    std::cout << "===== 神童射箭计算程序 Ver1.0 =====\n";
    std::cout << "===== 仅供学习娱乐！计算结果仅供参考 ==作者：bili猫盒喵===\n";
    std::cout << "输入 r(m) 和 alpha(度)，空格分隔，输入 000 退出程序\n";
    while (true) {
        std::cout << "\n请输入 r alpha: ";
        if (!std::getline(std::cin, line) || line == "000") break;
        std::istringstream iss(line);
        double r, alpha;
        if (!(iss >> r >> alpha)) {
            std::cout << "错误: 请输入两个数值 (r alpha)\n";
            continue;
        }
        double alpha_rad = deg2rad(alpha);

        double theta_low_orig = calculate_original_theta(r, alpha_rad, alpha, 20);
        double theta_high_orig = calculate_original_theta(r, alpha_rad, alpha, 60);
        double low_corr = calculate_correction_low(low_params, r, alpha);
        double high_corr = calculate_correction_high(high_params, r, alpha);

        double theta_low_corr = theta_low_orig + low_corr;
        double theta_high_corr = theta_high_orig + high_corr;
        // NaN handling
        if (std::isnan(theta_high_orig)) theta_high_corr = 90.0 - theta_low_corr;
        if (std::isnan(theta_low_orig))  theta_low_corr = 90.0 - theta_high_corr;

        theta_low_corr = constrain_theta(theta_low_corr);
        theta_high_corr = constrain_theta(theta_high_corr);

        std::cout << "\n计算结果:\n";
        std::cout << "r = " << r << " m, alpha = " << alpha << "°\n";
        std::cout << "----------------------------------------\n";
        std::cout << "低角度 theta: 原始=" << theta_low_orig
            << "°, 校正=" << theta_low_corr
            << "°, 校正量=" << (theta_low_corr - theta_low_orig) << "°\n";
        std::cout << "高角度 theta: 原始=" << theta_high_orig
            << "°, 校正=" << theta_high_corr
            << "°, 校正量=" << (theta_high_corr - theta_high_orig) << "°\n";
        std::cout << "----------------------------------------\n";

        bool warn = (std::fabs(theta_low_corr) > 90 || std::fabs(theta_high_corr) > 90);
        if (warn) {
            std::cout << "警告: 检测到角度值超出±90°范围!\n";
            if (std::fabs(theta_low_corr) > 90) {
                std::cout << "建议低角度替代值: " << (90.0 - theta_high_corr) << "°\n";
            }
            if (std::fabs(theta_high_corr) > 90) {
                std::cout << "建议高角度替代值: " << (90.0 - theta_low_corr) << "°\n";
            }
        }

        std::cout << "\n推荐使用值:\n";
        if (!warn) {
            std::cout << "  低角度 theta = " << theta_low_corr << "°\n";
            std::cout << "  高角度 theta = " << theta_high_corr - 2 << "°\n";
        }
        else {
            if (std::fabs(theta_low_corr) > 90)
                std::cout << "  低角度 theta = " << (90.0 - theta_high_corr) << "° (替代)\n";
            else
                std::cout << "  低角度 theta = " << theta_low_corr << "°\n";
            if (std::fabs(theta_high_corr) > 90)
                std::cout << "  高角度 theta = " << (90.0 - theta_low_corr) << "° (替代)\n";
            else
                std::cout << "  高角度 theta = " << theta_high_corr << "°\n";
        }
    }
    std::cout << "程序退出\n";
    return 0;
}
