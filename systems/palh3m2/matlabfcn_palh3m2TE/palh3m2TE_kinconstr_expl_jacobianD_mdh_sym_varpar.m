% Jacobian time derivative of explicit kinematic constraints of
% palh3m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% WD [12x4]
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)
% [ParkChoPlo1999] Park, FC and Choi, Jihyeon and Ploen, SR: Symbolic formulation of closed chain dynamics in independent coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function WD = palh3m2TE_kinconstr_expl_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_kinconstr_expl_jacobianD_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_kinconstr_expl_jacobianD_mdh_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_kinconstr_expl_jacobianD_mdh_sym_varpar: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 01:41:27
% EndTime: 2020-05-07 01:41:28
% DurationCPUTime: 1.05s
% Computational Cost: add. (4422->75), mult. (8188->144), div. (218->20), fcn. (10780->14), ass. (0->92)
t197 = sin(qJ(2));
t198 = sin(pkin(15));
t239 = cos(qJ(2));
t240 = cos(pkin(15));
t184 = t197 * t198 - t239 * t240;
t186 = t197 * t240 + t239 * t198;
t196 = sin(qJ(3));
t200 = cos(qJ(3));
t161 = t184 * t200 + t186 * t196;
t191 = pkin(17) + pkin(18);
t189 = sin(t191);
t190 = cos(t191);
t226 = t184 * t196;
t218 = t186 * t200 - t226;
t129 = t189 * t161 - t190 * t218;
t192 = sin(pkin(16));
t193 = cos(pkin(16));
t253 = t190 * t161 + t189 * t218;
t112 = t129 * t192 - t193 * t253;
t107 = 0.1e1 / t112 ^ 2;
t269 = t129 * t193 + t192 * t253;
t277 = t269 ^ 2 * t107;
t180 = t186 * qJD(2);
t179 = t184 * qJD(2);
t228 = t179 * t196;
t135 = -t228 - qJD(3) * t226 + (qJD(3) * t186 + t180) * t200;
t212 = -t161 * qJD(3) - t179 * t200 - t180 * t196;
t247 = t190 * t212;
t115 = t189 * t135 - t247;
t249 = t189 * t212;
t216 = t190 * t135 + t249;
t102 = t115 * t192 - t193 * t216;
t106 = 0.1e1 / t112;
t267 = t102 * t106;
t276 = t267 * t277;
t275 = t102 * t112;
t100 = 0.1e1 + t277;
t101 = t115 * t193 + t192 * t216;
t237 = t107 * t269;
t274 = 0.2e1 * (t106 * t112 + t237 * t269) * (t101 * t237 - t276) / t100 ^ 2;
t272 = t269 * t101;
t271 = 0.2e1 * t276;
t127 = 0.1e1 / t253 ^ 2;
t126 = 0.1e1 / t253;
t241 = t216 * t126;
t236 = t127 * t241;
t259 = 0.2e1 * t129;
t270 = t129 * t236 * t259;
t199 = sin(pkin(14));
t201 = cos(pkin(14));
t219 = -t184 * t201 + t186 * t199;
t154 = 0.1e1 / t219 ^ 2;
t167 = t184 * t199 + t186 * t201;
t265 = t154 * t167 ^ 2;
t125 = t129 ^ 2;
t123 = t125 * t127 + 0.1e1;
t234 = t127 * t129;
t268 = 0.2e1 * (-t126 * t253 - t129 * t234) / t123 ^ 2 * (t115 * t234 - t125 * t236);
t194 = sin(pkin(18));
t195 = cos(pkin(18));
t157 = t184 * t195 + t186 * t194;
t150 = 0.1e1 / t157 ^ 2;
t220 = -t184 * t194 + t186 * t195;
t258 = t150 * t220 ^ 2;
t262 = t216 * t253;
t147 = -t179 * t199 - t180 * t201;
t153 = 0.1e1 / t219;
t250 = t147 * t153;
t261 = t250 * t265;
t257 = (t179 * t201 - t180 * t199) * t167;
t149 = 0.1e1 / t157;
t214 = t179 * t194 - t180 * t195;
t256 = t214 * t149;
t255 = t256 * t258;
t251 = (-t179 * t195 - t180 * t194) * t220;
t181 = 0.1e1 / t184 ^ 2;
t183 = t186 ^ 2;
t171 = t183 * t181 + 0.1e1;
t227 = t180 / t184 * t181;
t229 = t179 * t181;
t235 = (-t183 * t227 - t186 * t229) / t171 ^ 2;
t169 = 0.1e1 / t171;
t144 = 0.1e1 + t265;
t141 = 0.1e1 + t258;
t136 = -t218 * qJD(3) - t180 * t200 + t228;
t121 = 0.1e1 / t123;
t120 = t189 * t136 + t247;
t117 = -t190 * t136 + t249;
t98 = 0.1e1 / t100;
t95 = t268 + (-t270 + t117 * t126 + (-t262 + (t115 - t120) * t129) * t127) * t121;
t94 = t268 + (-t270 + t241 + (t115 * t259 - t262) * t127) * t121;
t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, t274 + ((t117 * t193 + t120 * t192) * t106 + t271 + (t275 - (t117 * t192 - t120 * t193) * t269 - t272) * t107) * t98, t274 + (-t267 + t271 + (-0.2e1 * t272 + t275) * t107) * t98, 0; 0, 0, 0, 0; 0, 0.2e1 * (-t153 * t219 - t265) / t144 ^ 2 * (-t154 * t257 - t261) + (t250 - 0.2e1 * t261 + (-t147 * t219 - 0.2e1 * t257) * t154) / t144, 0, 0; 0, 0.2e1 * (t149 * t157 + t258) / t141 ^ 2 * (t150 * t251 + t255) + (t256 - 0.2e1 * t255 + (-t157 * t214 - 0.2e1 * t251) * t150) / t141, 0, 0; 0, t95, t94, 0; 0, -0.2e1 * t235 + 0.2e1 * (-t169 * t229 + (-t169 * t227 - t181 * t235) * t186) * t186, 0, 0; 0, t95, t94, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
WD = t1;
