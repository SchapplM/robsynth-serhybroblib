% Calculate matrix of centrifugal and coriolis load on the joints for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh2m1DE_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1DE_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'palh2m1DE_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:06
% EndTime: 2020-05-02 23:52:10
% DurationCPUTime: 0.84s
% Computational Cost: add. (1124->129), mult. (1974->213), div. (0->0), fcn. (1212->18), ass. (0->97)
t75 = sin(qJ(1));
t79 = cos(qJ(1));
t158 = -t75 ^ 2 - t79 ^ 2;
t73 = sin(qJ(3));
t74 = sin(qJ(2));
t131 = t74 * t73;
t77 = cos(qJ(3));
t62 = t77 * pkin(3) + pkin(2);
t78 = cos(qJ(2));
t101 = -pkin(3) * t131 + t62 * t78;
t93 = pkin(1) + t101;
t33 = pkin(4) + t93;
t51 = t79 * rSges(6,1) - t75 * rSges(6,2);
t53 = rSges(6,1) * t75 + t79 * rSges(6,2);
t72 = sin(qJ(4));
t76 = cos(qJ(4));
t11 = -t33 * t75 - t72 * t51 - t53 * t76;
t136 = rSges(6,2) * t72;
t138 = rSges(6,1) * t76;
t54 = -t136 + t138;
t126 = t33 + t54;
t49 = -t72 * rSges(6,1) - t76 * rSges(6,2);
t12 = t126 * t79 + t75 * t49;
t157 = m(6) * (t11 * t75 - t12 * t79) + m(5) * t158 * (rSges(5,1) + t93);
t71 = qJ(2) + qJ(3);
t65 = sin(t71);
t66 = cos(t71);
t156 = rSges(4,1) * t65 + rSges(4,2) * t66;
t151 = m(6) * pkin(3);
t110 = t49 * t151;
t128 = t78 * t77;
t39 = t128 - t131;
t104 = t39 * t110;
t155 = qJD(3) * t104;
t150 = -t66 / 0.2e1;
t152 = m(6) / 0.4e1;
t69 = qJ(4) + qJ(2);
t148 = cos(t69);
t70 = -qJ(4) + qJ(2);
t147 = cos(t70);
t146 = sin(t69);
t145 = sin(t70);
t144 = m(6) * rSges(6,2);
t143 = rSges(3,3) * m(3);
t142 = rSges(4,3) * m(4);
t141 = rSges(5,3) * m(5);
t140 = pkin(2) * t74;
t91 = rSges(4,1) * t66 - rSges(4,2) * t65 + t78 * pkin(2) + pkin(1);
t135 = (-t75 * rSges(4,3) + t91 * t79) * t79;
t132 = t49 * t79;
t129 = t78 * t73;
t125 = t156 * t75;
t124 = m(6) * qJD(2);
t123 = m(6) * qJD(4);
t122 = Icges(3,4) * t74;
t13 = -t79 * rSges(4,3) - t91 * t75;
t36 = pkin(3) * t129 + t74 * t62;
t80 = -0.2e1 * Icges(4,4) * t66 * t150 + (-Icges(4,4) * t65 + 0.2e1 * (-Icges(4,1) + Icges(4,2)) * t150) * t65;
t1 = m(3) * t158 * (rSges(3,1) * t74 + rSges(3,2) * t78) * (rSges(3,1) * t78 - rSges(3,2) * t74 + pkin(1)) - t74 * (-Icges(3,1) * t78 + t122) / 0.2e1 + t74 * (-Icges(3,2) * t78 - t122) / 0.2e1 + m(4) * (t13 * (t75 * t140 + t125) + (-t156 - t140) * t135) + t80 - (-0.2e1 * Icges(3,4) * t78 + (-Icges(3,1) + Icges(3,2)) * t74) * t78 / 0.2e1 + t157 * t36;
t118 = t1 * qJD(1);
t37 = (t74 * t77 + t129) * pkin(3);
t2 = m(4) * (t13 * t125 - t135 * t156) + t80 + t157 * t37;
t117 = t2 * qJD(1);
t94 = -t51 * t76 + t53 * t72;
t3 = t11 * t94 + t12 * (-t54 * t75 + t132);
t116 = t3 * qJD(1);
t105 = qJ(3) + t69;
t97 = sin(t105);
t106 = qJ(3) + t70;
t98 = sin(t106);
t81 = -t98 / 0.4e1 - t97 / 0.4e1;
t100 = cos(t106);
t99 = cos(t105);
t82 = t100 / 0.4e1 - t99 / 0.4e1;
t4 = ((t36 * t72 / 0.2e1 + t82 * pkin(3) + (t147 / 0.4e1 - t148 / 0.4e1) * pkin(2)) * rSges(6,2) + (-t36 * t76 / 0.2e1 + t81 * pkin(3) + (-t145 / 0.4e1 - t146 / 0.4e1) * pkin(2)) * rSges(6,1)) * m(6);
t115 = t4 * qJD(1);
t111 = pkin(3) * t144;
t85 = t99 * t111;
t86 = t100 * t111;
t112 = rSges(6,1) * t151;
t87 = t97 * t112;
t88 = t98 * t112;
t6 = -t65 * (rSges(4,2) * t142 - Icges(4,6)) + (rSges(4,1) * t142 + pkin(3) * t141 - Icges(4,5)) * t66 + t85 / 0.2e1 + t86 / 0.2e1 + t87 / 0.2e1 - t88 / 0.2e1;
t114 = t6 * qJD(3);
t90 = -t138 / 0.2e1 + t136 / 0.2e1;
t84 = t90 * t37;
t7 = (t84 + (t81 * rSges(6,1) + t82 * rSges(6,2)) * pkin(3)) * m(6);
t113 = t7 * qJD(1);
t107 = t54 * t123;
t102 = t88 / 0.4e1 - t86 / 0.4e1 + t87 / 0.4e1 + t85 / 0.4e1;
t34 = ((m(4) * rSges(4,1) + pkin(3) * (m(5) + m(6))) * t73 + m(4) * t77 * rSges(4,2)) * pkin(2);
t9 = (-t39 / 0.2e1 - t131 / 0.2e1 + t128 / 0.2e1) * t110;
t83 = t34 * qJD(2) + t9 * qJD(4);
t26 = t34 * qJD(3);
t8 = m(6) * t84 + t102;
t5 = t90 * t36 * m(6) + t102 + (-t147 * t144 / 0.4e1 + (rSges(6,2) * t148 + (t145 + t146) * rSges(6,1)) * t152) * pkin(2);
t10 = [0.4e1 * (t11 * (-t33 * t79 + t94) + t12 * (-t126 * t75 + t132)) * t152 * qJD(1) + t1 * qJD(2) + t2 * qJD(3) + t3 * t123, t118 + (((t141 + t142) * pkin(2) + rSges(3,1) * t143 - Icges(3,5)) * t78 - t74 * (rSges(3,2) * t143 - Icges(3,6)) + t6) * qJD(2) + t114 + t5 * qJD(4) + ((t147 / 0.2e1 + t148 / 0.2e1) * rSges(6,2) + (-t145 / 0.2e1 + t146 / 0.2e1) * rSges(6,1)) * pkin(2) * t124, t6 * qJD(2) + t8 * qJD(4) + t114 + t117, t5 * qJD(2) + t8 * qJD(3) + (t49 * qJD(4) * t33 + t116) * m(6); -t4 * qJD(4) - t118, -t26, -t26 - t83, -t9 * qJD(3) + t36 * t107 - t115; -t7 * qJD(4) - t117, t83, 0, t9 * qJD(2) + t37 * t107 - t113; -m(6) * t116 + t4 * qJD(2) + t7 * qJD(3), -t101 * t49 * t124 + t115 - t155, -qJD(2) * t104 + t113 - t155, 0;];
Cq = t10;
