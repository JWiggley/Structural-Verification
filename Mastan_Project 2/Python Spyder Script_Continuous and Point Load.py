import numpy as np
import matplotlib.pyplot as plt

def calculate_beam_reactions(length, uniform_load, point_load, uniform_end, point_loc):
    total_uniform = uniform_load * uniform_end
    uniform_center = uniform_end / 2.0
    
    reaction_a = (total_uniform * (length - uniform_center) + point_load * (length - point_loc)) / length
    reaction_b = (total_uniform + point_load) - reaction_a
    
    return reaction_a, reaction_b

def analyze_span_profiles(length, E, I, wu, Pu, rect_end, point_loc):
    points = np.linspace(0, length, 1000)
    shear = np.zeros_like(points)
    moment = np.zeros_like(points)
    deflection = np.zeros_like(points)
    
    Rua, Rub = calculate_beam_reactions(length, wu, Pu, rect_end, point_loc)
    
    w_service = wu / 1.2
    P_service = Pu / 1.6
    R_service, _ = calculate_beam_reactions(length, w_service, P_service, rect_end, point_loc)
    
    for i, x in enumerate(points):
        w_dist = x if x < rect_end else rect_end
        p_active = 1.0 if x > point_loc else 0.0
        shear[i] = Rua - (wu * w_dist) - (Pu * p_active)
        
        w_moment = (wu * x**2 / 2.0) if x <= rect_end else (wu * rect_end * (x - rect_end / 2.0))
        p_moment = Pu * (x - point_loc) if x > point_loc else 0.0
        moment[i] = (Rua * x) - w_moment - p_moment
        
        mac_w1 = (x**4) if x > 0 else 0
        mac_w2 = ((x - rect_end)**4) if x > rect_end else 0
        mac_P = ((x - point_loc)**3) if x > point_loc else 0
        
        ei_y = ((R_service / 6.0) * x**3 - (w_service / 24.0) * mac_w1 
                + (w_service / 24.0) * mac_w2 - (P_service / 6.0) * mac_P - 22500.0 * x)
        deflection[i] = ei_y / (E * I)
        
    return points, shear, moment, deflection

def generate_design_plots(x, V, M, y):
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    
    ax1.plot(x, V, color='crimson', lw=2)
    ax1.fill_between(x, V, color='crimson', alpha=0.15)
    ax1.axhline(0, color='black', lw=0.8, ls='--')
    ax1.set_ylabel('Shear V (kips)')
    ax1.set_title('Beam Design Profile (W10x15)')
    ax1.grid(True, ls=':')
    
    ax2.plot(x, M, color='navy', lw=2)
    ax2.fill_between(x, M, color='navy', alpha=0.15)
    ax2.axhline(0, color='black', lw=0.8, ls='--')
    ax2.set_ylabel('Moment M (kip-in)')
    ax2.grid(True, ls=':')
    
    ax3.plot(x, y, color='darkgreen', lw=2)
    ax3.fill_between(x, y, color='darkgreen', alpha=0.15)
    ax3.axhline(0, color='black', lw=0.8, ls='--')
    ax3.set_xlabel('Beam Position x (inches)')
    ax3.set_ylabel('Deflection y (in)')
    ax3.grid(True, ls=':')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    span_length = 120.0  
    modulus = 29000.0    
    inertia = 68.9       
    
    factored_uniform = 0.50  
    factored_point = 16.0    
    uniform_span = 60.0      
    point_position = 90.0   
    
    x, shear, moment, deflection = analyze_span_profiles(
        span_length, modulus, inertia, 
        factored_uniform, factored_point, 
        uniform_span, point_position
    )
    
    generate_design_plots(x, shear, moment, deflection)
    
    Rua, Rub = calculate_beam_reactions(span_length, factored_uniform, factored_point, uniform_span, point_position)
    print(f"Left Reaction: {Rua:.2f} kips")
    print(f"Right Reaction: {Rub:.2f} kips")
    print(f"Max Moment: {np.max(moment):.2f} kip-in")
    print(f"Max Deflection: {np.min(deflection):.3f} inches")
