import numpy

class freevortex:
    def calc_hub_angles(r_m, r_hub, alpha_2_rad, alpha_3_rad):
        alpha_2_hub_rad = numpy.arctan((r_m/r_hub) *numpy.tan(alpha_2_rad))
        alpha_2_hub_deg = numpy.rad2deg(alpha_2_hub_rad)
    
        alpha_3_hub_rad = numpy.arctan((r_m/r_hub) *numpy.tan(alpha_3_rad))
        alpha_3_hub_deg = numpy.rad2deg(alpha_3_hub_rad)
    
        beta_2_hub_rad = numpy.arctan((r_m/r_hub) *numpy.tan(alpha_2_rad) - (r/r_m)*(flow_coeff_2)**(-1))
        beta_2_hub_deg = numpy.rad2deg(beta_2_hub_rad)
    
        beta_3_hub_rad = numpy.arctan((r_m/r_hub) *numpy.tan(alpha_3_rad) - (r/r_m)*(flow_coeff_3)**(-1))
        beta_3_hub_deg = numpy.rad2deg(beta_3_hub_rad)  
    
        return
    
    def calc_tip_angles():
        alpha_2_tip_rad = numpy.arctan((r_m/r_tip) *numpy.tan(alpha_2_rad))
        alpha_2_tip_deg = numpy.rad2deg(alpha_2_hub_tip)
    
        alpha_3_tip_rad = numpy.arctan((r_m/r_tip) *numpy.tan(alpha_3_rad))
        alpha_3_tip_deg = numpy.rad2deg(alpha_3_tip_rad)
    
        beta_2_tip_rad = numpy.arctan((r_m/r_tip) *numpy.tan(alpha_2_rad) - (r_tip/r_m)*(flow_coeff_2)**(-1))
        beta_2_tip_deg = numpy.rad2deg(beta_2_tip_rad)
    
        beta_3_tip_rad = numpy.arctan((r_m/r_tip) *numpy.tan(alpha_3_rad) - (r_tip/r_m)*(flow_coeff_3)**(-1))
        beta_3_tip_deg = numpy.rad2deg(beta_3_tip_rad)      
        return
    
    def check_reaction_hub():
    
        R2 = (Ca_3 * numpy.tan(beta_3_hub) - Ca_2*numpy.tan(beta_2_hub))/(2*U_hub)
    
        return
    
    
    def check_reaction_tip():
    
        R2 = (Ca_3 * numpy.tan(beta_3_tip) - Ca_2*numpy.tan(beta_2_tip))/(2*U_tip)
    
        return