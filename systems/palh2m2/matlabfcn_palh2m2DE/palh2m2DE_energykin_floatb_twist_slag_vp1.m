% Calculate kinetic energy for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = palh2m2DE_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2DE_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:22
% EndTime: 2020-05-03 01:06:24
% DurationCPUTime: 2.16s
% Computational Cost: add. (565->255), mult. (928->338), div. (0->0), fcn. (686->8), ass. (0->126)
t211 = Icges(2,4) - Icges(4,5) - Icges(6,5);
t210 = Icges(2,1) + Icges(4,1) + Icges(6,1);
t209 = Icges(4,4) + Icges(6,4) + Icges(2,5);
t208 = Icges(2,2) + Icges(4,3) + Icges(6,3);
t207 = Icges(2,6) - Icges(4,6) - Icges(6,6);
t137 = sin(qJ(1));
t206 = t211 * t137;
t141 = cos(qJ(1));
t205 = t211 * t141;
t126 = V_base(6) + qJD(1);
t140 = cos(qJ(2));
t131 = t140 * pkin(4);
t132 = pkin(1) + pkin(2);
t115 = t131 + t132;
t139 = cos(qJ(3));
t189 = t139 * pkin(5);
t152 = t115 + t189;
t204 = t126 * (pkin(3) + t152);
t133 = rSges(4,1) + pkin(1);
t203 = -t133 - t131;
t202 = -t141 * t208 - t206;
t201 = t137 * t208 - t205;
t198 = t137 * t210 + t205;
t197 = t141 * t210 - t206;
t135 = sin(qJ(3));
t159 = t139 * rSges(5,1) - t135 * rSges(5,2);
t196 = (t115 + t159) * t126;
t136 = sin(qJ(2));
t160 = t140 * rSges(3,1) - t136 * rSges(3,2);
t150 = pkin(1) + t160;
t195 = -t141 * rSges(3,3) + t150 * t137;
t194 = t137 * rSges(3,3) + t150 * t141;
t117 = V_base(6) * t131;
t127 = rSges(6,1) + t132;
t193 = (t127 + t189) * V_base(6) + (rSges(6,1) + t152) * qJD(1) + t117;
t191 = pkin(4) * t136;
t190 = pkin(5) * t135;
t188 = pkin(4) * qJD(2);
t187 = pkin(5) * qJD(3);
t134 = sin(qJ(4));
t138 = cos(qJ(4));
t76 = t141 * t134 + t137 * t138;
t186 = Icges(7,4) * t76;
t184 = t137 * rSges(4,3);
t183 = t137 * rSges(5,3);
t180 = Icges(3,4) * t136;
t179 = Icges(3,4) * t140;
t178 = Icges(5,4) * t135;
t177 = Icges(5,4) * t139;
t119 = V_base(4) * t131;
t174 = V_base(4) * t189 + t119;
t173 = V_base(5) * t191 + V_base(1);
t172 = t140 * t188 + V_base(3);
t169 = t136 * t188;
t168 = t135 * t187;
t166 = V_base(5) * t190 + t173;
t165 = t139 * t187 + t172;
t163 = -t190 - t191;
t162 = -t131 - t189;
t161 = t133 * V_base(6) + t117;
t114 = qJD(2) * t137 + V_base(4);
t158 = Icges(3,1) * t140 - t180;
t157 = Icges(5,1) * t139 - t178;
t156 = -Icges(3,2) * t136 + t179;
t155 = -Icges(5,2) * t135 + t177;
t154 = Icges(3,5) * t140 - Icges(3,6) * t136;
t153 = Icges(5,5) * t139 - Icges(5,6) * t135;
t151 = t132 + t159;
t148 = -t168 - t169;
t111 = -qJD(3) * t141 + V_base(5);
t113 = qJD(3) * t137 + V_base(4);
t147 = (-Icges(5,3) * t141 + t153 * t137) * t111 + (Icges(5,3) * t137 + t153 * t141) * t113 + (Icges(5,5) * t135 + Icges(5,6) * t139) * t126;
t112 = -qJD(2) * t141 + V_base(5);
t146 = (-Icges(3,3) * t141 + t154 * t137) * t112 + (Icges(3,3) * t137 + t154 * t141) * t114 + (Icges(3,5) * t136 + Icges(3,6) * t140) * t126;
t145 = -t114 * t191 + V_base(2);
t144 = t126 * rSges(6,3) + t148;
t65 = -Icges(5,6) * t141 + t155 * t137;
t66 = Icges(5,6) * t137 + t155 * t141;
t69 = -Icges(5,5) * t141 + t157 * t137;
t70 = Icges(5,5) * t137 + t157 * t141;
t91 = Icges(5,2) * t139 + t178;
t99 = Icges(5,1) * t135 + t177;
t143 = (-t135 * t66 + t139 * t70) * t113 + (-t135 * t65 + t139 * t69) * t111 + (-t135 * t91 + t139 * t99) * t126;
t102 = Icges(3,1) * t136 + t179;
t67 = -Icges(3,6) * t141 + t156 * t137;
t68 = Icges(3,6) * t137 + t156 * t141;
t71 = -Icges(3,5) * t141 + t158 * t137;
t72 = Icges(3,5) * t137 + t158 * t141;
t94 = Icges(3,2) * t140 + t180;
t142 = (-t136 * t68 + t140 * t72) * t114 + (-t136 * t67 + t140 * t71) * t112 + (t102 * t140 - t136 * t94) * t126;
t125 = pkin(3) + t132;
t122 = qJD(4) + t126;
t110 = t141 * rSges(2,1) - t137 * rSges(2,2);
t109 = t137 * rSges(2,1) + t141 * rSges(2,2);
t108 = t136 * rSges(3,1) + t140 * rSges(3,2);
t107 = t135 * rSges(5,1) + t139 * rSges(5,2);
t106 = V_base(5) * rSges(7,1) - rSges(7,2) * V_base(4);
t105 = rSges(7,1) * V_base(4) + V_base(5) * rSges(7,2);
t80 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t79 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t78 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t77 = -t137 * t134 + t141 * t138;
t73 = Icges(7,4) * t77;
t60 = V_base(5) * rSges(2,3) - t126 * t109 + V_base(1);
t59 = -V_base(4) * rSges(2,3) + t126 * t110 + V_base(2);
t58 = V_base(4) * t109 - V_base(5) * t110 + V_base(3);
t57 = Icges(7,1) * t77 - t186;
t56 = Icges(7,1) * t76 + t73;
t55 = -Icges(7,2) * t76 + t73;
t54 = Icges(7,2) * t77 + t186;
t51 = (-rSges(4,3) * V_base(4) + t203 * V_base(5)) * t141 + (-rSges(4,3) * V_base(5) + t133 * V_base(4) + t119) * t137 + t172;
t50 = t112 * t108 - t195 * t126 + V_base(1);
t49 = -t114 * t108 + t194 * t126 + V_base(2);
t48 = V_base(5) * rSges(4,2) + (qJD(1) * t203 - t161) * t137 + (t126 * rSges(4,3) - t169) * t141 + t173;
t47 = (-t203 * t141 + t184) * qJD(1) + t161 * t141 + V_base(6) * t184 - V_base(4) * rSges(4,2) + t145;
t46 = t160 * qJD(2) - t194 * V_base(5) + t195 * V_base(4) + V_base(3);
t45 = -rSges(6,3) * V_base(4) * t141 + (t127 * V_base(4) + t174) * t137 + ((-t127 + t162) * t141 - t137 * rSges(6,3)) * V_base(5) + t165;
t44 = -t113 * t107 + t126 * t183 + t196 * t141 + t145;
t43 = t111 * t107 + (t126 * rSges(5,3) - t169) * t141 - t196 * t137 + t173;
t42 = t159 * qJD(3) + t119 * t137 + (-t141 * rSges(5,3) + t151 * t137) * V_base(4) + ((-t151 - t131) * t141 - t183) * V_base(5) + t172;
t41 = V_base(5) * rSges(6,2) - t193 * t137 + t144 * t141 + t166;
t40 = V_base(2) + t193 * t141 + (-rSges(6,2) + t163) * V_base(4) + t144 * t137;
t39 = V_base(2) + (t122 * (rSges(7,1) * t138 - rSges(7,2) * t134) + t204) * t141 + ((-t134 * rSges(7,1) - t138 * rSges(7,2)) * t122 + t148) * t137 + (-rSges(7,3) + t163) * V_base(4);
t38 = -t141 * t169 - t141 * t168 + V_base(5) * rSges(7,3) - t122 * ((rSges(7,1) * t137 + t141 * rSges(7,2)) * t138 + t134 * (t141 * rSges(7,1) - rSges(7,2) * t137)) - t137 * t204 + t166;
t37 = (t105 * t138 + t106 * t134 + t125 * V_base(4) + t174) * t137 + t165 + (t105 * t134 - t106 * t138 + (-t125 + t162) * V_base(5)) * t141;
t1 = m(1) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + t114 * (t146 * t137 + t142 * t141) / 0.2e1 + t112 * (t142 * t137 - t146 * t141) / 0.2e1 + t113 * (t147 * t137 + t143 * t141) / 0.2e1 + t111 * (t143 * t137 - t147 * t141) / 0.2e1 + m(2) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(3) * (t46 ^ 2 + t49 ^ 2 + t50 ^ 2) / 0.2e1 + m(4) * (t47 ^ 2 + t48 ^ 2 + t51 ^ 2) / 0.2e1 + m(7) * (t37 ^ 2 + t38 ^ 2 + t39 ^ 2) / 0.2e1 + m(5) * (t42 ^ 2 + t43 ^ 2 + t44 ^ 2) / 0.2e1 + m(6) * (t40 ^ 2 + t41 ^ 2 + t45 ^ 2) / 0.2e1 + ((t136 * t72 + t140 * t68) * t114 + (t136 * t71 + t140 * t67) * t112 + (t135 * t70 + t139 * t66) * t113 + (t135 * t69 + t139 * t65) * t111 + (t136 * t102 + t135 * t99 + t139 * t91 + t140 * t94 + Icges(4,2) + Icges(6,2) + Icges(2,3)) * t126) * t126 / 0.2e1 + ((t202 * t137 + t198 * t141 - t76 * t54 + t77 * t56 + Icges(1,4)) * V_base(5) + (t201 * t137 + t197 * t141 - t76 * t55 + t77 * t57 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t198 * t137 - t202 * t141 + t77 * t54 + t76 * t56 + Icges(1,2)) * V_base(5) + (t197 * t137 - t201 * t141 + t77 * t55 + t76 * t57 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t126 * (t137 * t209 + t141 * t207) + V_base(4) * t126 * (-t207 * t137 + t209 * t141) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(7,5) * t76 + Icges(7,6) * t77) * V_base(5) + (Icges(7,5) * t77 - Icges(7,6) * t76) * V_base(4) + Icges(7,3) * t122 / 0.2e1) * t122;
T = t1;
