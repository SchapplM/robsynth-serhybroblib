% Jacobian of implicit kinematic constraints of
% picker2Dm2IC
% projection from active to passive joints coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% B21 [(no of passive joints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 09:21
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function B21 = picker2Dm2IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: qJ has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_projection_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 09:20:36
% EndTime: 2020-05-11 09:20:36
% DurationCPUTime: 0.25s
% Computational Cost: add. (627->100), mult. (194->94), div. (48->10), fcn. (219->58), ass. (0->71)
t166 = qJ(1) + qJ(2);
t158 = sin(t166);
t159 = cos(t166);
t168 = sin(qJ(7));
t171 = cos(qJ(7));
t203 = (t158 * t168 + t159 * t171) * pkin(3);
t196 = qJ(6) - qJ(3);
t186 = -qJ(9) + t196;
t187 = qJ(4) - t196;
t188 = -qJ(4) - t196;
t120 = 0.1e1 / ((-cos(t188) + cos(t187)) * pkin(6) + (cos(qJ(4) + t186) - cos(qJ(4) - t186)) * pkin(2));
t202 = pkin(2) * t120;
t201 = cos(qJ(1)) * pkin(1);
t173 = 0.1e1 / pkin(5);
t190 = (pkin(8) + qJ(5));
t153 = t190 - t196;
t165 = qJ(1) + qJ(8);
t176 = t153 - t165;
t154 = t190 + t196;
t177 = t154 - t165;
t200 = 0.1e1 / (pkin(6) * (cos(t177) - cos(t176)) + (-cos(qJ(9) - t177) + cos(qJ(9) + t176)) * pkin(2)) * t173;
t175 = 0.1e1 / pkin(3);
t199 = t120 * t175;
t197 = qJ(2) + qJ(4);
t162 = qJ(1) + t197;
t151 = sin(t162);
t152 = cos(t162);
t180 = t151 * t159 - t158 * t152;
t121 = 0.1e1 / t180;
t174 = 0.1e1 / pkin(4);
t198 = t121 * t174;
t150 = -t188 + t166;
t138 = -qJ(9) + t150;
t149 = t187 + t166;
t139 = qJ(9) + t149;
t195 = sin(t139) - sin(t138);
t194 = cos(t139) - cos(t138);
t193 = sin(t149) - sin(t150);
t192 = cos(t149) - cos(t150);
t191 = -cos(0.2e1 * t165) + cos((2 * t190));
t155 = qJ(9) + t162;
t156 = -qJ(9) + t162;
t124 = sin(t156) - sin(t155);
t125 = cos(t156) - cos(t155);
t169 = sin(qJ(4));
t164 = 0.1e1 / t169;
t189 = pkin(1) * t164 * t174;
t185 = -qJ(1) + t190;
t184 = -qJ(4) + t185;
t183 = qJ(4) + t185;
t182 = -qJ(8) + t184;
t181 = -qJ(8) + t183;
t179 = t151 * t168 + t152 * t171;
t170 = sin(qJ(2));
t167 = sin(qJ(8));
t163 = sin(qJ(1)) * pkin(1);
t161 = -qJ(9) + t190;
t160 = qJ(9) + t190;
t157 = sin(t197);
t148 = qJ(9) + t153;
t147 = -qJ(9) + t154;
t140 = -0.1e1 / sin(qJ(8) - t185);
t129 = -pkin(5) * sin(t165) + t163;
t128 = t201 - pkin(5) * cos(t165);
t127 = cos(t161) - cos(t160);
t126 = sin(t161) - sin(t160);
t123 = pkin(3) * t159 + pkin(4) * t152 + t201;
t122 = pkin(3) * t158 + pkin(4) * t151 + t163;
t119 = t124 * t171 - t125 * t168;
t118 = t179 * t121;
t1 = [(-t157 * pkin(1) - t169 * pkin(3)) * t164 * t175, -t118; ((-t195 * t122 - t194 * t123) * t199 + ((sin(t148) - sin(t147)) * t129 + (cos(t148) - cos(t147)) * t128) * t200) * pkin(2), (-t194 * t168 + t195 * t171) * t202; (t170 * pkin(3) + t157 * pkin(4)) * t175 * t189, (t179 * pkin(4) + t203) * t198; t167 * t140, 0; ((t124 * t122 + t125 * t123) * t199 + (-t126 * t129 - t128 * t127) * t200) * pkin(2), -t119 * t202; (-t191 * pkin(5) + (-cos(qJ(8) + 0.2e1 * qJ(1)) + cos((2 * pkin(8)) + (2 * qJ(5)) + qJ(8))) * pkin(1)) * t173 / t191, 0; -(pkin(3) * (cos(t184) - cos(t183)) + (cos(-qJ(2) + t182) - cos(qJ(2) + t181)) * pkin(5)) * pkin(1) * t175 * t173 / (cos(t182) - cos(t181)), -t118; -t170 * t189, (t180 * pkin(4) - t203) * t198; (pkin(5) * t167 - pkin(1) * sin(t185)) * t173 * t140, 0; ((pkin(2) * t126 + (-sin(t154) + sin(t153)) * pkin(6)) * t129 + t128 * (pkin(2) * t127 + (-cos(t154) + cos(t153)) * pkin(6))) * t200 + (-(pkin(2) * t124 + t193 * pkin(6)) * t122 - (pkin(2) * t125 + t192 * pkin(6)) * t123) * t199, (pkin(2) * t119 + (-t192 * t168 + t193 * t171) * pkin(6)) * t120;];
B21 = t1;
