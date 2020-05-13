% Calculate kinetic energy for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m1DE_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1DE_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:51:58
% EndTime: 2020-05-02 23:52:00
% DurationCPUTime: 1.96s
% Computational Cost: add. (614->240), mult. (901->322), div. (0->0), fcn. (699->10), ass. (0->122)
t207 = Icges(2,4) + Icges(5,5);
t206 = Icges(2,1) + Icges(5,1);
t205 = -Icges(5,4) + Icges(2,5);
t204 = Icges(2,2) + Icges(5,3);
t203 = Icges(2,6) - Icges(5,6);
t133 = cos(qJ(3));
t118 = t133 * pkin(3) + pkin(2);
t130 = sin(qJ(2));
t134 = cos(qJ(2));
t129 = sin(qJ(3));
t184 = pkin(3) * t129;
t202 = t118 * t134 - t130 * t184;
t135 = cos(qJ(1));
t201 = t207 * t135;
t131 = sin(qJ(1));
t200 = t207 * t131;
t120 = V_base(6) + qJD(1);
t137 = pkin(1) + rSges(5,1);
t156 = pkin(1) + t202;
t199 = (rSges(5,1) + t156) * qJD(1) + (t137 + t202) * V_base(6);
t196 = t120 * (pkin(4) + t156);
t194 = -rSges(6,3) - pkin(5) - pkin(6);
t193 = t204 * t135 + t200;
t192 = -t204 * t131 + t201;
t191 = t206 * t131 + t201;
t190 = t206 * t135 - t200;
t171 = t134 * t129;
t172 = t130 * t118;
t75 = pkin(3) * t171 + t172;
t121 = qJD(2) * t135;
t114 = V_base(5) + t121;
t159 = -t134 * rSges(3,1) + t130 * rSges(3,2);
t149 = pkin(1) - t159;
t189 = t135 * rSges(3,3) + t131 * t149;
t188 = -t131 * rSges(3,3) + t135 * t149;
t185 = pkin(2) * t130;
t183 = t134 * pkin(2);
t128 = sin(qJ(4));
t132 = cos(qJ(4));
t85 = t135 * t128 + t131 * t132;
t182 = Icges(6,4) * t85;
t177 = Icges(3,4) * t130;
t176 = Icges(3,4) * t134;
t127 = qJ(2) + qJ(3);
t122 = sin(t127);
t175 = Icges(4,4) * t122;
t123 = cos(t127);
t174 = Icges(4,4) * t123;
t170 = V_base(5) * pkin(5) + V_base(1);
t163 = t134 * V_base(4);
t169 = t163 * t184 + V_base(2);
t168 = pkin(3) * (t130 * t133 + t171) * qJD(3);
t162 = t120 * (pkin(1) + t183);
t115 = -qJD(2) * t131 + V_base(4);
t158 = rSges(4,1) * t123 - rSges(4,2) * t122;
t157 = -V_base(4) * pkin(5) + V_base(2);
t155 = Icges(3,1) * t134 - t177;
t154 = Icges(4,1) * t123 - t175;
t153 = -Icges(3,2) * t130 + t176;
t152 = -Icges(4,2) * t122 + t174;
t151 = Icges(3,5) * t134 - Icges(3,6) * t130;
t150 = Icges(4,5) * t123 - Icges(4,6) * t122;
t148 = -t75 * qJD(2) - t168;
t82 = qJD(3) * t135 + t114;
t83 = V_base(4) + (-qJD(2) - qJD(3)) * t131;
t147 = (-Icges(4,5) * t122 - Icges(4,6) * t123) * t120 + (Icges(4,3) * t135 + t131 * t150) * t82 + (-Icges(4,3) * t131 + t135 * t150) * t83;
t146 = (Icges(3,3) * t135 + t131 * t151) * t114 + (-Icges(3,3) * t131 + t135 * t151) * t115 + (-Icges(3,5) * t130 - Icges(3,6) * t134) * t120;
t145 = t131 * V_base(4) - V_base(5) * t135;
t143 = -qJD(2) + t145;
t144 = -pkin(3) * (qJD(3) - t143) * t123 - qJD(2) * t183 + V_base(3);
t142 = -rSges(5,3) * t120 + t148;
t62 = Icges(4,6) * t135 + t131 * t152;
t63 = -Icges(4,6) * t131 + t135 * t152;
t64 = Icges(4,5) * t135 + t131 * t154;
t65 = -Icges(4,5) * t131 + t135 * t154;
t79 = -Icges(4,2) * t123 - t175;
t80 = -Icges(4,1) * t122 - t174;
t141 = (-t122 * t63 + t123 * t65) * t83 + (-t122 * t62 + t123 * t64) * t82 + (-t122 * t79 + t123 * t80) * t120;
t104 = -Icges(3,1) * t130 - t176;
t70 = Icges(3,6) * t135 + t131 * t153;
t71 = -Icges(3,6) * t131 + t135 * t153;
t72 = Icges(3,5) * t135 + t131 * t155;
t73 = -Icges(3,5) * t131 + t135 * t155;
t99 = -Icges(3,2) * t134 - t177;
t140 = (-t130 * t71 + t134 * t73) * t115 + (-t130 * t70 + t134 * t72) * t114 + (t104 * t134 - t130 * t99) * t120;
t139 = pkin(1) + pkin(4);
t136 = pkin(5) - rSges(5,2);
t117 = qJD(4) + t120;
t116 = pkin(2) * t163;
t111 = t135 * rSges(2,1) - t131 * rSges(2,2);
t110 = t131 * rSges(2,1) + t135 * rSges(2,2);
t109 = -t130 * rSges(3,1) - t134 * rSges(3,2);
t108 = V_base(5) * rSges(6,1) - rSges(6,2) * V_base(4);
t107 = rSges(6,1) * V_base(4) + V_base(5) * rSges(6,2);
t91 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t90 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t89 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t86 = -t131 * t128 + t135 * t132;
t81 = -t122 * rSges(4,1) - t123 * rSges(4,2);
t77 = Icges(6,4) * t86;
t67 = -t131 * rSges(4,3) + t135 * t158;
t66 = t135 * rSges(4,3) + t131 * t158;
t59 = V_base(5) * rSges(2,3) - t120 * t110 + t170;
t58 = t120 * t111 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t57 = V_base(4) * t110 - V_base(5) * t111 + V_base(3);
t56 = Icges(6,1) * t86 - t182;
t55 = Icges(6,1) * t85 + t77;
t54 = -Icges(6,2) * t85 + t77;
t53 = Icges(6,2) * t86 + t182;
t50 = t114 * t109 - t189 * t120 + t170;
t49 = -t115 * t109 + t188 * t120 + t157;
t48 = t159 * qJD(2) - t188 * V_base(5) + t189 * V_base(4) + V_base(3);
t47 = (V_base(4) * rSges(5,3) + (-t137 - t183) * V_base(5)) * t135 + (V_base(5) * rSges(5,3) + t137 * V_base(4) + t116) * t131 + t144;
t46 = t115 * t185 + t120 * t67 + t135 * t162 - t83 * t81 + t157;
t45 = -t114 * t185 - t120 * t66 - t162 * t131 + t82 * t81 + t170;
t44 = pkin(1) * t145 + t143 * t183 + t83 * t66 - t82 * t67 + V_base(3);
t43 = (t107 * t128 - t108 * t132 + (-t139 - t183) * V_base(5)) * t135 + (t107 * t132 + t108 * t128 + t139 * V_base(4) + t116) * t131 + t144;
t42 = V_base(1) - t199 * t131 + t142 * t135 + (t136 - t75) * V_base(5);
t41 = t199 * t135 + (-t136 + t172) * V_base(4) + t142 * t131 + t169;
t40 = (t117 * (rSges(6,1) * t132 - rSges(6,2) * t128) + t196) * t135 + (t117 * (-t128 * rSges(6,1) - t132 * rSges(6,2)) + t148) * t131 + (t172 + t194) * V_base(4) + t169;
t39 = -t75 * t121 - t135 * t168 + V_base(1) - t117 * ((rSges(6,1) * t131 + t135 * rSges(6,2)) * t132 + t128 * (t135 * rSges(6,1) - t131 * rSges(6,2))) - t131 * t196 + (-t75 - t194) * V_base(5);
t1 = m(1) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(2) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + m(3) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) / 0.2e1 + t115 * (-t146 * t131 + t140 * t135) / 0.2e1 + t114 * (t140 * t131 + t146 * t135) / 0.2e1 + m(4) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + t83 * (-t147 * t131 + t135 * t141) / 0.2e1 + t82 * (t131 * t141 + t135 * t147) / 0.2e1 + m(5) * (t41 ^ 2 + t42 ^ 2 + t47 ^ 2) / 0.2e1 + m(6) * (t39 ^ 2 + t40 ^ 2 + t43 ^ 2) / 0.2e1 + ((-t130 * t73 - t134 * t71) * t115 + (-t130 * t72 - t134 * t70) * t114 + (-t122 * t65 - t123 * t63) * t83 + (-t122 * t64 - t123 * t62) * t82 + (-t130 * t104 - t122 * t80 - t123 * t79 - t134 * t99 + Icges(5,2) + Icges(2,3)) * t120) * t120 / 0.2e1 + ((-t193 * t131 + t191 * t135 - t85 * t53 + t86 * t55 + Icges(1,4)) * V_base(5) + (-t192 * t131 + t190 * t135 - t85 * t54 + t86 * t56 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t191 * t131 + t193 * t135 + t86 * t53 + t85 * t55 + Icges(1,2)) * V_base(5) + (t190 * t131 + t192 * t135 + t86 * t54 + t85 * t56 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t120 * (t205 * t131 + t203 * t135) + V_base(4) * t120 * (-t203 * t131 + t205 * t135) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(6,5) * t85 + Icges(6,6) * t86) * V_base(5) + (Icges(6,5) * t86 - Icges(6,6) * t85) * V_base(4) + Icges(6,3) * t117 / 0.2e1) * t117;
T = t1;
