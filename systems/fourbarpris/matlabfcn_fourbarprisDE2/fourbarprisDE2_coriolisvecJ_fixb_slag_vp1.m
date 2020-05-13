% Calculate vector of centrifugal and Coriolis load on the joints for
% fourbarprisDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% m [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [1x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbarprisDE2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_coriolisvecJ_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE2_coriolisvecJ_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_coriolisvecJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE2_coriolisvecJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisDE2_coriolisvecJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisDE2_coriolisvecJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:45:03
% EndTime: 2020-05-07 09:45:11
% DurationCPUTime: 2.23s
% Computational Cost: add. (4449->114), mult. (4438->245), div. (314->16), fcn. (228->2), ass. (0->116)
t74 = qJ(1) + pkin(3);
t163 = -pkin(2) - t74;
t53 = pkin(1) - t163;
t42 = 1 / t53;
t162 = -pkin(2) + t74;
t55 = pkin(1) + t162;
t46 = 1 / t55;
t173 = t42 * t46;
t56 = pkin(1) + t163;
t189 = 1 / t56;
t54 = pkin(1) - t162;
t190 = 1 / t54;
t196 = t189 * t190;
t146 = t173 * t196;
t64 = t74 ^ 2;
t170 = (0.1e1 / t64 ^ 2 / pkin(1) ^ 2);
t65 = 0.1e1 / t74;
t130 = (t65 * t146 * t170);
t144 = (t46 * t189 * t170);
t45 = 1 / (t54 ^ 2);
t47 = 1 / (t55 ^ 2);
t49 = 1 / (t56 ^ 2);
t201 = 4 * t130 + (-t45 * t144 + (t189 * t47 - t46 * t49) * t190 * t170) * t42;
t182 = 2 * qJ(1);
t85 = qJ(1) ^ 2;
t185 = 3 * t85;
t88 = (pkin(3) ^ 2);
t184 = 3 * t88;
t89 = pkin(2) ^ 2;
t90 = pkin(1) ^ 2;
t63 = (-t89 - t90);
t87 = pkin(3) * t88;
t95 = t88 ^ 2;
t26 = 2 * t87 * rSges(3,3) + t63 * t184 + 3 * t95;
t27 = 3 * pkin(3) * t63 + rSges(3,3) * t184 + 7 * t87;
t75 = pkin(3) * rSges(3,3);
t33 = 2 * t75 + 8 * t88 + t63;
t78 = rSges(3,3) / 0.2e1;
t59 = t78 + 0.9e1 / 0.2e1 * pkin(3);
t84 = qJ(1) * t85;
t98 = t85 ^ 2;
t199 = (t27 * t182) + (t33 * t185) + 0.4e1 * t59 * t84 + t26 + (5 * t98);
t172 = t55 * t56;
t159 = t54 * t172;
t143 = t53 * t159;
t93 = sqrt(-t143);
t19 = 0.1e1 / t93;
t198 = -t19 / 0.2e1;
t76 = pkin(1) + pkin(2);
t77 = pkin(1) - pkin(2);
t169 = t76 ^ 2 * t77 ^ 2;
t139 = -2 * pkin(3) * qJ(1) - t85 - t88;
t168 = (t90 - t89);
t31 = t139 + t168;
t176 = rSges(3,1) * t31;
t179 = t95 / 0.2e1;
t5 = -t93 * t176 / 0.2e1 + t59 * t98 + (t33 * t84) + (t27 * t85) + (t63 * t87) + ((t98 + t26) * qJ(1)) + (t179 - t169 / 0.2e1) * rSges(3,3) + (t179 + t169 / 0.2e1) * pkin(3);
t4 = t5 * qJD(1);
t197 = qJD(1) * t65;
t37 = (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3);
t195 = t37 * t196;
t52 = t75 + 2 * t88;
t60 = 0.5e1 / 0.2e1 * pkin(3) + t78;
t193 = (t182 * t60 + t185 + t52) * t93;
t138 = t190 * t144;
t125 = qJD(1) * t138;
t43 = 0.1e1 / t53 ^ 2;
t192 = t201 * qJD(1) + t43 * t125;
t91 = 0.1e1 / pkin(1);
t188 = t91 ^ 2;
t39 = t74 * qJD(1);
t187 = -0.2e1 * t39;
t186 = -0.2e1 * t74;
t181 = t19 / 0.2e1;
t180 = t31 / 0.2e1;
t12 = -t159 + (t172 + (t55 - t56) * t54) * t53;
t11 = t12 * qJD(1);
t151 = t19 / t143 * t180;
t66 = 0.1e1 / t74 ^ 2;
t174 = t31 * t66;
t161 = t19 * t174;
t81 = qJD(1) ^ 2;
t178 = (-t81 * t161 + (t11 * t151 + t187 * t19) * t197) * t65;
t177 = (t65 * t12 * t151 + (t186 * t65 - t174) * t19) * t197;
t175 = t19 * t31;
t67 = t65 / t64;
t171 = t67 * t81;
t167 = qJD(1) * t31;
t158 = t19 * t167;
t157 = t66 * t167;
t156 = rSges(2,1) * t198;
t155 = rSges(2,2) * t181;
t150 = t90 / 0.2e1 - t89 / 0.2e1;
t17 = t84 + t60 * t85 + (t52 * qJ(1)) + t87 / 0.2e1 + t150 * pkin(3) + (t88 / 0.2e1 - t150) * rSges(3,3);
t154 = t17 * t181;
t153 = -t175 / 0.2e1;
t152 = -t175 / 0.4e1;
t141 = qJD(1) * t153;
t140 = t158 / 0.2e1;
t30 = -t139 + t168;
t136 = 0.2e1 * rSges(3,1) * (t180 - t30 / 0.2e1);
t135 = t66 * t141;
t134 = t66 * t140;
t131 = t42 * t138;
t10 = t17 * t93 + (t30 * t176) / 0.2e1;
t29 = t31 ^ 2;
t122 = qJD(1) * t11 * t29 * t43 * t45 * t47 * t49;
t121 = t42 * t125;
t113 = t178 / 0.2e1 + t153 * t171;
t107 = (qJD(1) * t193 + t11 * t154 + t136 * t39) * t131;
t16 = -rSges(2,1) * t93 - rSges(2,2) * t30;
t15 = -rSges(2,1) * t30 + rSges(2,2) * t93;
t9 = t10 * qJD(1);
t7 = rSges(2,2) * t187 + t11 * t156;
t6 = rSges(2,1) * t187 + t11 * t155;
t1 = [m(2) * ((-(-t113 * t16 + t135 * t7) * t16 + (t113 * t15 + t134 * t6) * t15) * t188 * t161 - (t15 * t6 + t16 * t7) * t29 * t121 / 0.2e1 + (t15 ^ 2 + t16 ^ 2) * (t31 * t39 * t121 + t29 * t81 * t130 - t122 * t170 / 0.4e1)) / 0.2e1 + m(2) * (-t16 * ((rSges(2,2) * t186 + t12 * t156) * t135 + (-t177 / 0.2e1 + t67 * t140) * t16) + ((rSges(2,1) * t186 + t12 * t155) * t134 + (t177 / 0.2e1 + t67 * t141) * t15) * t15) * t157 * t188 * t198 + 0.2e1 * ((t42 * t47) + t43 * t46) * t64 * t81 * t195 + (-t107 * t9 + (-qJD(1) * t107 + t192 * t9) * t10 + (t192 * t5 - 0.2e1 * ((t11 * t152 + t39 * t93) * rSges(3,1) + t199 * qJD(1)) * t131) * t4 + (t9 * (t12 * t154 + t74 * t136 + t193) + t4 * ((t12 * t152 + t74 * t93) * rSges(3,1) + t199)) * t121 - (t4 ^ 2 + t9 ^ 2) * (t43 * t138 + t201) / 0.2e1) * m(3) + (Icges(3,2) + Icges(2,3)) * (t175 * t178 - t66 * t122 / 0.2e1 - t158 * t177 + (0.2e1 * t39 * t157 + t29 * t171) * t146) + (-0.4e1 * t74 * t195 - 0.2e1 * (t189 * t45 + t190 * t49) * t37 * t64) * t81 * t173;];
tauc = t1(:);
