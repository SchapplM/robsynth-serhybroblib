% Calculate homogenous joint transformation matrices for
% palh3m1DE1
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
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh3m1DE1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 10:23:45
% EndTime: 2020-04-18 10:23:47
% DurationCPUTime: 1.74s
% Computational Cost: add. (26288->80), mult. (39422->124), div. (1920->6), fcn. (24950->42), ass. (0->111)
t282 = -pkin(6) - pkin(2);
t281 = -pkin(6) + pkin(2);
t280 = -pkin(8) - pkin(10);
t279 = -pkin(8) + pkin(10);
t245 = pkin(6) ^ 2;
t250 = pkin(2) ^ 2;
t231 = sin(qJ(2));
t233 = sin(pkin(16));
t237 = cos(qJ(2));
t239 = cos(pkin(16));
t218 = t231 * t233 - t237 * t239;
t271 = pkin(5) * t218;
t258 = pkin(1) * t271;
t217 = -0.2e1 * t258;
t252 = pkin(1) ^ 2;
t261 = pkin(5) ^ 2 + t252;
t256 = t217 + t261;
t213 = -t245 + t250 + t256;
t215 = pkin(1) - t271;
t262 = t217 + t252;
t211 = sqrt(-((pkin(5) - t281) * (pkin(5) + t281) + t262) * ((pkin(5) - t282) * (pkin(5) + t282) + t262));
t219 = t231 * t239 + t237 * t233;
t268 = t211 * t219;
t205 = -pkin(5) * t268 + t215 * t213;
t278 = -t205 / 0.2e1;
t206 = pkin(5) * t219 * t213 + t215 * t211;
t277 = t206 / 0.2e1;
t223 = sin(pkin(17));
t276 = t223 / 0.2e1;
t225 = sin(pkin(19));
t275 = t225 / 0.2e1;
t230 = sin(qJ(3));
t274 = t230 / 0.2e1;
t273 = cos(pkin(15)) / 0.2e1;
t236 = cos(qJ(3));
t214 = 0.1e1 / t256;
t251 = 0.1e1 / pkin(2);
t266 = t214 * t251;
t200 = (t206 * t274 + t236 * t278) * t266;
t201 = (t205 * t274 + t236 * t277) * t266;
t222 = pkin(18) + pkin(19);
t220 = sin(t222);
t221 = cos(t222);
t193 = -t221 * t200 - t220 * t201;
t272 = pkin(4) * t193;
t259 = pkin(3) * t272;
t191 = -0.2e1 * t259;
t248 = pkin(4) ^ 2;
t263 = t191 + t248;
t185 = sqrt(-((pkin(3) - t279) * (pkin(3) + t279) + t263) * ((pkin(3) - t280) * (pkin(3) + t280) + t263));
t192 = t220 * t200 - t221 * t201;
t270 = t185 * t192;
t260 = pkin(3) ^ 2 + t248;
t257 = t191 + t260;
t188 = 0.1e1 / t257;
t242 = 0.1e1 / pkin(10);
t269 = t188 * t242;
t246 = 0.1e1 / pkin(6);
t267 = t214 * t246;
t244 = 0.1e1 / pkin(8);
t265 = t242 * t244;
t264 = t246 * t251;
t255 = -t250 + t261;
t243 = pkin(8) ^ 2;
t254 = -t243 + t260;
t253 = t188 * t244 / 0.2e1;
t241 = pkin(10) ^ 2;
t238 = cos(qJ(1));
t235 = cos(qJ(4));
t234 = sin(pkin(15));
t232 = sin(qJ(1));
t229 = sin(qJ(4));
t228 = cos(pkin(18));
t227 = cos(pkin(19));
t226 = sin(pkin(18));
t224 = cos(pkin(17));
t216 = pkin(1) * t218 - pkin(5);
t212 = t217 + t245 + t255;
t210 = atan2(t211 * t264 / 0.2e1, -(t245 - t255 + 0.2e1 * t258) * t264 / 0.2e1);
t209 = cos(t210);
t208 = sin(t210);
t207 = pkin(1) * t219 * t212 - t216 * t211;
t204 = -pkin(1) * t268 - t216 * t212;
t203 = t225 * t208 - t227 * t209;
t202 = -t227 * t208 - t225 * t209;
t199 = atan2((t207 * t273 - t204 * t234 / 0.2e1) * t267, (t204 * t273 + t207 * t234 / 0.2e1) * t267);
t198 = atan2((t205 * t275 + t227 * t277) * t266, (t206 * t275 + t227 * t278) * t266);
t197 = cos(t199);
t196 = sin(t199);
t195 = cos(t198);
t194 = sin(t198);
t190 = -pkin(3) * t193 + pkin(4);
t189 = -pkin(3) + t272;
t187 = t191 + t241 + t254;
t186 = -t241 + t243 + t257;
t184 = atan2(t185 * t265 / 0.2e1, -(t241 - t254 + 0.2e1 * t259) * t265 / 0.2e1);
t183 = cos(t184);
t182 = sin(t184);
t181 = pkin(3) * t192 * t187 + t190 * t185;
t180 = -pkin(3) * t270 + t190 * t187;
t179 = -t223 * t182 + t224 * t183;
t178 = t224 * t182 + t223 * t183;
t177 = atan2((pkin(4) * t192 * t186 - t189 * t185) * t253, (-pkin(4) * t270 - t189 * t186) * t253);
t176 = cos(t177);
t175 = sin(t177);
t174 = atan2((t181 * t224 / 0.2e1 + t180 * t276) * t269, (-t180 * t224 / 0.2e1 + t181 * t276) * t269);
t173 = cos(t174);
t172 = sin(t174);
t171 = -t226 * t175 - t228 * t176;
t170 = -t228 * t175 + t226 * t176;
t1 = [t238, -t232, 0, 0; t232, t238, 0, 0; 0, 0, 1, pkin(12); 0, 0, 0, 1; t237, -t231, 0, pkin(13); 0, 0, -1, 0; t231, t237, 0, 0; 0, 0, 0, 1; -t236, t230, 0, pkin(1); -t230, -t236, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t173, -t172, 0, pkin(4); t172, t173, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t235, -t229, 0, pkin(9); 0, 0, -1, -pkin(11); t229, t235, 0, 0; 0, 0, 0, 1; t197, -t196, 0, -pkin(7); 0, 0, -1, 0; t196, t197, 0, pkin(14); 0, 0, 0, 1; t195, -t194, 0, pkin(1); t194, t195, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t171, -t170, 0, t228 * pkin(3); t170, t171, 0, -t226 * pkin(3); 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(6); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -1, 0, 0, pkin(8); 0, -1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t203, -t202, 0, t227 * pkin(2); t202, t203, 0, t225 * pkin(2); 0, 0, 1, 0; 0, 0, 0, 1; t179, -t178, 0, t224 * pkin(10); t178, t179, 0, t223 * pkin(10); 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,12);             % numerisch
else,                         T_mdh = sym('xx', [4,4,12]); end % symbolisch

for i = 1:12
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
