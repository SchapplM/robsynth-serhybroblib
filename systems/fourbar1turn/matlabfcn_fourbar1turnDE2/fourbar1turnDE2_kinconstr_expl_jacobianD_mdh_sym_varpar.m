% Jacobian time derivative of explicit kinematic constraints of
% fourbar1turnDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% WD [6x2]
%
% Sources:
% [NakamuraGho1989] Nakamura, Yoshihiko and Ghodoussi, Modjtaba: Dynamics computation of closed-link robot mechanisms with nonredundant and redundant actuators (1989)
% [ParkChoPlo1999] Park, FC and Choi, Jihyeon and Ploen, SR: Symbolic formulation of closed chain dynamics in independent coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function WD = fourbar1turnDE2_kinconstr_expl_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_kinconstr_expl_jacobianD_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_kinconstr_expl_jacobianD_mdh_sym_varpar: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_kinconstr_expl_jacobianD_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_expl_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:33:16
% EndTime: 2020-04-12 19:33:18
% DurationCPUTime: 1.27s
% Computational Cost: add. (4790->117), mult. (6672->262), div. (188->18), fcn. (1744->4), ass. (0->125)
t148 = cos(qJ(2));
t226 = pkin(2) * t148;
t198 = pkin(1) * t226;
t145 = -0.2e1 * t198;
t150 = pkin(4) ^ 2;
t152 = pkin(3) ^ 2;
t154 = pkin(2) ^ 2;
t155 = pkin(1) ^ 2;
t205 = t154 + t155;
t189 = -t152 + t205;
t138 = t145 + t150 + t189;
t143 = pkin(1) - t226;
t147 = sin(qJ(2));
t206 = t145 + t155;
t234 = -pkin(3) - pkin(4);
t131 = (pkin(2) - t234) * (pkin(2) + t234) + t206;
t233 = -pkin(3) + pkin(4);
t132 = (pkin(2) - t233) * (pkin(2) + t233) + t206;
t215 = t131 * t132;
t156 = sqrt(-t215);
t210 = t147 * t156;
t113 = -pkin(2) * t210 + t138 * t143;
t108 = 0.1e1 / t113 ^ 2;
t227 = pkin(2) * t147;
t129 = t138 * t227;
t114 = t143 * t156 + t129;
t110 = t114 ^ 2;
t102 = t108 * t110 + 0.1e1;
t220 = t108 * t114;
t107 = 0.1e1 / t113;
t177 = pkin(1) * pkin(2) * (-t131 - t132);
t121 = t147 * t177;
t120 = qJD(2) * t121;
t125 = 0.1e1 / t156;
t219 = t125 * t120;
t191 = t147 * t219;
t200 = 0.2e1 * t143 * pkin(1);
t209 = t148 * t156;
t90 = (-t191 + (-t209 + (t138 + t200) * t147) * qJD(2)) * pkin(2);
t222 = t107 * t108 * t90;
t146 = t147 ^ 2;
t228 = pkin(1) * t154;
t194 = t146 * t228;
t179 = qJD(2) * t194;
t203 = qJD(2) * t156;
t188 = t147 * t203;
t204 = qJD(2) * t148;
t208 = (t138 * t204 + t188) * pkin(2);
t217 = t125 * t143;
t91 = t120 * t217 + 0.2e1 * t179 + t208;
t238 = (-t110 * t222 + t91 * t220) / t102 ^ 2;
t218 = t125 * t121;
t142 = t145 + t205;
t137 = -t150 + t152 + t142;
t144 = pkin(1) * t148 - pkin(2);
t112 = -pkin(1) * t210 - t137 * t144;
t104 = 0.1e1 / t112;
t190 = 0.2e1 * t198;
t136 = t150 - t189 + t190;
t133 = 0.1e1 / t136;
t139 = 0.1e1 / t142;
t105 = 0.1e1 / t112 ^ 2;
t134 = 0.1e1 / t136 ^ 2;
t140 = 0.1e1 / t142 ^ 2;
t236 = 0.2e1 * t125;
t235 = -0.2e1 * t144;
t232 = -t139 / 0.2e1;
t231 = t139 / 0.2e1;
t230 = pkin(1) * t137;
t229 = pkin(1) * t140;
t225 = pkin(3) * t142;
t224 = pkin(4) * t142;
t202 = pkin(2) * t235;
t89 = (-t191 + (-t209 + (t137 + t202) * t147) * qJD(2)) * pkin(1);
t223 = t104 * t105 * t89;
t130 = t147 * t230;
t115 = -t144 * t156 + t130;
t221 = t105 * t115;
t216 = t125 * t144;
t214 = t132 * t134;
t213 = t139 * t148;
t212 = t140 * t147;
t211 = t146 * t155;
t207 = pkin(1) * t188 + t204 * t230;
t201 = 0.2e1 * t224;
t199 = pkin(1) * t227;
t197 = 0.2e1 * t140;
t98 = 0.1e1 / t102;
t196 = t98 * t224;
t195 = t104 * t225;
t193 = t139 * t219;
t192 = t154 * t211;
t187 = qJD(2) * t211;
t186 = t134 * t199;
t185 = t140 * t199;
t184 = qJD(2) * t199;
t151 = 0.1e1 / pkin(4);
t182 = t151 * t196;
t180 = 0.6e1 * t147 * t204;
t178 = pkin(2) * t187;
t176 = 0.4e1 / t215 * t120 * t218;
t175 = -0.2e1 * t185;
t153 = 0.1e1 / pkin(3);
t93 = t130 + (-t209 + (t202 - t218) * t147) * pkin(1);
t86 = (t112 * t175 + t139 * t93) * t153;
t174 = 0.2e1 * t115 * t86 * t225;
t172 = t139 * t140 * t154 * t187;
t171 = -t176 / 0.4e1;
t170 = t176 / 0.4e1;
t168 = 0.8e1 * t172;
t167 = -0.2e1 * pkin(4) * t98 * t184;
t116 = (t148 * t177 - 0.4e1 * t192) * qJD(2);
t166 = (t147 * t171 + (-t147 * t116 / 0.2e1 - t148 * t120) * t236) * t139;
t135 = t133 * t134;
t119 = -t131 * t214 + 0.1e1;
t111 = t115 ^ 2;
t103 = t105 * t111 + 0.1e1;
t96 = t121 * t217 + 0.2e1 * t194 + (t138 * t148 + t210) * pkin(2);
t95 = -t121 * t216 + 0.2e1 * pkin(2) * t211 + (t137 * t148 + t210) * pkin(1);
t94 = t129 + (-t209 + (t200 - t218) * t147) * pkin(2);
t92 = -t120 * t216 + 0.2e1 * t178 + t207;
t88 = (t115 * t175 + t139 * t95) * t153;
t87 = (-t114 * t185 + t96 * t231) * t151;
t85 = (t113 * t185 + t94 * t232) * t151;
t1 = [0, 0; 0, 0; 0, (t105 * t174 - 0.2e1 * t88 * t195) / t103 ^ 2 * (-t111 * t223 + t92 * t221) + (0.2e1 * (t88 * t104 - t86 * t221) * pkin(3) * t184 + ((t155 * pkin(2) * t180 - t116 * t216 + t144 * t171) * t139 + t115 * t168 + ((-0.2e1 * t92 * t140 * pkin(2) + t193) * t147 + ((t209 + (-t137 + t218) * t147) * t139 + (-t115 * t148 - t147 * t95) * pkin(2) * t197) * qJD(2)) * pkin(1)) * t153 * t195 + t174 * t223 + (-t88 * t89 - ((0.4e1 * t178 + t207) * t139 + t112 * t168 + (t166 + (-0.2e1 * t89 * t212 + (t213 * t235 + (-t112 * t148 - t147 * t93) * t197) * qJD(2)) * pkin(2)) * pkin(1)) * t153 * t115 - t86 * t92) * t105 * t225) / t103; 0, 0.2e1 * (-((0.4e1 * t179 + t208) * t232 - 0.4e1 * t113 * t172 + (-t166 / 0.2e1 + (t90 * t212 + (-t143 * t213 + (t113 * t148 + t147 * t94) * t140) * qJD(2)) * pkin(1)) * pkin(2)) * t182 + t85 * t167) * t220 + 0.2e1 * (-t85 * t91 + t87 * t90) * t108 * t196 + 0.2e1 * (t108 * t238 + t98 * t222) * t114 * t85 * t201 + 0.2e1 * (-((t116 * t217 + t143 * t170 + t180 * t228) * t231 + 0.4e1 * t114 * t172 + ((t193 / 0.2e1 - t91 * t229) * t147 + ((t209 + (-t138 + t218) * t147) * t231 + (-t148 * t114 - t147 * t96) * t229) * qJD(2)) * pkin(2)) * t182 + (t201 * t238 + t167) * t87) * t107; 0, (-0.2e1 * t133 * t218 - 0.4e1 * t156 * t186) * (-t214 + (-0.2e1 * t132 * t135 - t134) * t131) / t119 ^ 2 * t184 + (t133 * t170 + (t134 * t190 + 0.8e1 * t135 * t192) * t203 + (t116 * t133 / 0.2e1 + 0.2e1 * t120 * t186) * t236) / t119; 0, 0;];
WD = t1;
