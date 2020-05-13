% Calculate homogenous joint transformation matrices for
% palh3m1TE
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
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh3m1TE_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [19x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:17:20
% EndTime: 2020-04-17 15:17:21
% DurationCPUTime: 0.98s
% Computational Cost: add. (13153->80), mult. (19726->126), div. (960->6), fcn. (12462->24), ass. (0->107)
t275 = -pkin(6) - pkin(2);
t274 = -pkin(6) + pkin(2);
t273 = -pkin(8) - pkin(10);
t272 = -pkin(8) + pkin(10);
t212 = sin(pkin(17));
t271 = -t212 / 0.2e1;
t270 = t212 / 0.2e1;
t213 = cos(pkin(17));
t269 = -t213 / 0.2e1;
t268 = t213 / 0.2e1;
t214 = sin(pkin(19));
t267 = t214 / 0.2e1;
t216 = cos(pkin(19));
t266 = -t216 / 0.2e1;
t265 = t216 / 0.2e1;
t217 = cos(pkin(18));
t264 = -t217 / 0.2e1;
t219 = sin(qJ(3));
t263 = t219 / 0.2e1;
t262 = cos(pkin(15)) / 0.2e1;
t234 = pkin(6) ^ 2;
t239 = pkin(2) ^ 2;
t220 = sin(qJ(2));
t222 = sin(pkin(16));
t226 = cos(qJ(2));
t228 = cos(pkin(16));
t207 = t220 * t222 - t226 * t228;
t260 = pkin(5) * t207;
t246 = pkin(1) * t260;
t206 = -0.2e1 * t246;
t241 = pkin(1) ^ 2;
t249 = pkin(5) ^ 2 + t241;
t244 = t206 + t249;
t202 = -t234 + t239 + t244;
t204 = pkin(1) - t260;
t250 = t206 + t241;
t199 = sqrt(-((pkin(5) - t274) * (pkin(5) + t274) + t250) * ((pkin(5) - t275) * (pkin(5) + t275) + t250));
t208 = t220 * t228 + t226 * t222;
t256 = t199 * t208;
t194 = -pkin(5) * t256 + t204 * t202;
t195 = pkin(5) * t208 * t202 + t204 * t199;
t225 = cos(qJ(3));
t203 = 0.1e1 / t244;
t240 = 0.1e1 / pkin(2);
t254 = t203 * t240;
t190 = (-t194 * t225 / 0.2e1 + t195 * t263) * t254;
t191 = (t195 * t225 / 0.2e1 + t194 * t263) * t254;
t211 = pkin(18) + pkin(19);
t209 = sin(t211);
t210 = cos(t211);
t186 = -t210 * t190 - t209 * t191;
t261 = pkin(4) * t186;
t247 = pkin(3) * t261;
t184 = -0.2e1 * t247;
t237 = pkin(4) ^ 2;
t251 = t184 + t237;
t177 = sqrt(-((pkin(3) - t272) * (pkin(3) + t272) + t251) * ((pkin(3) - t273) * (pkin(3) + t273) + t251));
t185 = t209 * t190 - t210 * t191;
t259 = t177 * t185;
t248 = pkin(3) ^ 2 + t237;
t245 = t184 + t248;
t181 = 0.1e1 / t245;
t231 = 0.1e1 / pkin(10);
t258 = t181 * t231;
t233 = 0.1e1 / pkin(8);
t257 = t181 * t233;
t235 = 0.1e1 / pkin(6);
t255 = t203 * t235;
t253 = t231 * t233;
t252 = t235 * t240;
t243 = -t239 + t249;
t232 = pkin(8) ^ 2;
t242 = -t232 + t248;
t230 = pkin(10) ^ 2;
t227 = cos(qJ(1));
t224 = cos(qJ(4));
t223 = sin(pkin(15));
t221 = sin(qJ(1));
t218 = sin(qJ(4));
t215 = sin(pkin(18));
t205 = pkin(1) * t207 - pkin(5);
t201 = t206 + t234 + t243;
t200 = t234 - t243 + 0.2e1 * t246;
t198 = (t199 * t267 + t200 * t265) * t252;
t197 = (t199 * t266 + t200 * t267) * t252;
t196 = pkin(1) * t208 * t201 - t205 * t199;
t193 = -pkin(1) * t256 - t205 * t201;
t192 = (t196 * t262 - t193 * t223 / 0.2e1) * t255;
t189 = (t193 * t262 + t196 * t223 / 0.2e1) * t255;
t188 = (t194 * t267 + t195 * t265) * t254;
t187 = (t194 * t266 + t195 * t267) * t254;
t183 = -pkin(3) * t186 + pkin(4);
t182 = -pkin(3) + t261;
t180 = t184 + t230 + t242;
t179 = -t230 + t232 + t245;
t178 = t230 - t242 + 0.2e1 * t247;
t176 = (t177 * t271 + t178 * t269) * t253;
t175 = (t177 * t268 + t178 * t271) * t253;
t174 = pkin(3) * t185 * t180 + t183 * t177;
t173 = pkin(4) * t185 * t179 - t182 * t177;
t172 = -pkin(3) * t259 + t183 * t180;
t171 = -pkin(4) * t259 - t182 * t179;
t170 = (t171 * t264 - t215 * t173 / 0.2e1) * t257;
t169 = (t215 * t171 / 0.2e1 + t173 * t264) * t257;
t168 = (t172 * t270 + t174 * t268) * t258;
t167 = (t172 * t269 + t174 * t270) * t258;
t1 = [t227, -t221, 0, 0; t221, t227, 0, 0; 0, 0, 1, pkin(12); 0, 0, 0, 1; t226, -t220, 0, pkin(13); 0, 0, -1, 0; t220, t226, 0, 0; 0, 0, 0, 1; -t225, t219, 0, pkin(1); -t219, -t225, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t167, -t168, 0, pkin(4); t168, t167, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t224, -t218, 0, pkin(9); 0, 0, -1, -pkin(11); t218, t224, 0, 0; 0, 0, 0, 1; t189, -t192, 0, -pkin(7); 0, 0, -1, 0; t192, t189, 0, pkin(14); 0, 0, 0, 1; t187, -t188, 0, pkin(1); t188, t187, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t170, -t169, 0, t217 * pkin(3); t169, t170, 0, -t215 * pkin(3); 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(6); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -1, 0, 0, pkin(8); 0, -1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t198, -t197, 0, t216 * pkin(2); t197, t198, 0, t214 * pkin(2); 0, 0, 1, 0; 0, 0, 0, 1; t176, -t175, 0, t213 * pkin(10); t175, t176, 0, t212 * pkin(10); 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,12);             % numerisch
else,                         T_mdh = sym('xx', [4,4,12]); end % symbolisch

for i = 1:12
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
