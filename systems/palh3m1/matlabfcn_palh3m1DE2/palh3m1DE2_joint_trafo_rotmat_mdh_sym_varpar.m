% Calculate homogenous joint transformation matrices for
% palh3m1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% T_mdh [4x4x12]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh3m1DE2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:33:06
% EndTime: 2020-04-19 19:33:08
% DurationCPUTime: 1.78s
% Computational Cost: add. (26288->80), mult. (39422->124), div. (1920->6), fcn. (24950->42), ass. (0->111)
t263 = -pkin(6) - pkin(2);
t262 = -pkin(6) + pkin(2);
t261 = -pkin(8) - pkin(10);
t260 = -pkin(8) + pkin(10);
t226 = pkin(6) ^ 2;
t231 = pkin(2) ^ 2;
t212 = sin(qJ(2));
t214 = sin(pkin(16));
t218 = cos(qJ(2));
t220 = cos(pkin(16));
t199 = t212 * t214 - t218 * t220;
t252 = pkin(5) * t199;
t239 = pkin(1) * t252;
t198 = -0.2e1 * t239;
t233 = pkin(1) ^ 2;
t242 = pkin(5) ^ 2 + t233;
t237 = t198 + t242;
t194 = -t226 + t231 + t237;
t196 = pkin(1) - t252;
t243 = t198 + t233;
t192 = sqrt(-((pkin(5) - t262) * (pkin(5) + t262) + t243) * ((pkin(5) - t263) * (pkin(5) + t263) + t243));
t200 = t212 * t220 + t218 * t214;
t249 = t192 * t200;
t186 = -pkin(5) * t249 + t196 * t194;
t259 = -t186 / 0.2e1;
t187 = pkin(5) * t200 * t194 + t196 * t192;
t258 = t187 / 0.2e1;
t204 = sin(pkin(17));
t257 = t204 / 0.2e1;
t206 = sin(pkin(19));
t256 = t206 / 0.2e1;
t211 = sin(qJ(3));
t255 = t211 / 0.2e1;
t254 = cos(pkin(15)) / 0.2e1;
t217 = cos(qJ(3));
t195 = 0.1e1 / t237;
t232 = 0.1e1 / pkin(2);
t247 = t195 * t232;
t181 = (t187 * t255 + t217 * t259) * t247;
t182 = (t186 * t255 + t217 * t258) * t247;
t203 = pkin(18) + pkin(19);
t201 = sin(t203);
t202 = cos(t203);
t174 = -t202 * t181 - t201 * t182;
t253 = pkin(4) * t174;
t240 = pkin(3) * t253;
t172 = -0.2e1 * t240;
t229 = pkin(4) ^ 2;
t244 = t172 + t229;
t166 = sqrt(-((pkin(3) - t260) * (pkin(3) + t260) + t244) * ((pkin(3) - t261) * (pkin(3) + t261) + t244));
t173 = t201 * t181 - t202 * t182;
t251 = t166 * t173;
t241 = pkin(3) ^ 2 + t229;
t238 = t172 + t241;
t169 = 0.1e1 / t238;
t223 = 0.1e1 / pkin(10);
t250 = t169 * t223;
t227 = 0.1e1 / pkin(6);
t248 = t195 * t227;
t225 = 0.1e1 / pkin(8);
t246 = t223 * t225;
t245 = t227 * t232;
t236 = -t231 + t242;
t224 = pkin(8) ^ 2;
t235 = -t224 + t241;
t234 = t169 * t225 / 0.2e1;
t222 = pkin(10) ^ 2;
t219 = cos(qJ(1));
t216 = cos(qJ(4));
t215 = sin(pkin(15));
t213 = sin(qJ(1));
t210 = sin(qJ(4));
t209 = cos(pkin(18));
t208 = cos(pkin(19));
t207 = sin(pkin(18));
t205 = cos(pkin(17));
t197 = pkin(1) * t199 - pkin(5);
t193 = t198 + t226 + t236;
t191 = atan2(t192 * t245 / 0.2e1, -(t226 - t236 + 0.2e1 * t239) * t245 / 0.2e1);
t190 = cos(t191);
t189 = sin(t191);
t188 = pkin(1) * t200 * t193 - t197 * t192;
t185 = -pkin(1) * t249 - t197 * t193;
t184 = t206 * t189 - t208 * t190;
t183 = -t208 * t189 - t206 * t190;
t180 = atan2((t188 * t254 - t185 * t215 / 0.2e1) * t248, (t185 * t254 + t188 * t215 / 0.2e1) * t248);
t179 = atan2((t186 * t256 + t208 * t258) * t247, (t187 * t256 + t208 * t259) * t247);
t178 = cos(t180);
t177 = sin(t180);
t176 = cos(t179);
t175 = sin(t179);
t171 = -pkin(3) * t174 + pkin(4);
t170 = -pkin(3) + t253;
t168 = t172 + t222 + t235;
t167 = -t222 + t224 + t238;
t165 = atan2(t166 * t246 / 0.2e1, -(t222 - t235 + 0.2e1 * t240) * t246 / 0.2e1);
t164 = cos(t165);
t163 = sin(t165);
t162 = pkin(3) * t173 * t168 + t171 * t166;
t161 = -pkin(3) * t251 + t171 * t168;
t160 = -t204 * t163 + t205 * t164;
t159 = t205 * t163 + t204 * t164;
t158 = atan2((pkin(4) * t173 * t167 - t170 * t166) * t234, (-pkin(4) * t251 - t170 * t167) * t234);
t157 = cos(t158);
t156 = sin(t158);
t155 = atan2((t162 * t205 / 0.2e1 + t161 * t257) * t250, (-t161 * t205 / 0.2e1 + t162 * t257) * t250);
t154 = cos(t155);
t153 = sin(t155);
t152 = -t207 * t156 - t209 * t157;
t151 = -t209 * t156 + t207 * t157;
t1 = [t219, -t213, 0, 0; t213, t219, 0, 0; 0, 0, 1, pkin(12); 0, 0, 0, 1; t218, -t212, 0, pkin(13); 0, 0, -1, 0; t212, t218, 0, 0; 0, 0, 0, 1; -t217, t211, 0, pkin(1); -t211, -t217, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t154, -t153, 0, pkin(4); t153, t154, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t216, -t210, 0, pkin(9); 0, 0, -1, -pkin(11); t210, t216, 0, 0; 0, 0, 0, 1; t178, -t177, 0, -pkin(7); 0, 0, -1, 0; t177, t178, 0, pkin(14); 0, 0, 0, 1; t176, -t175, 0, pkin(1); t175, t176, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t152, -t151, 0, t209 * pkin(3); t151, t152, 0, -t207 * pkin(3); 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(6); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -1, 0, 0, pkin(8); 0, -1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t184, -t183, 0, t208 * pkin(2); t183, t184, 0, t206 * pkin(2); 0, 0, 1, 0; 0, 0, 0, 1; t160, -t159, 0, t205 * pkin(10); t159, t160, 0, t204 * pkin(10); 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,12);             % numerisch
else,                         T_mdh = sym('xx', [4,4,12]); end % symbolisch

for i = 1:12
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
