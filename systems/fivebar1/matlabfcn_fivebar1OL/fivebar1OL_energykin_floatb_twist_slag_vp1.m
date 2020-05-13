% Calculate kinetic energy for
% fivebar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fivebar1OL_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fivebar1OL_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:12:55
% EndTime: 2020-04-27 06:12:56
% DurationCPUTime: 0.60s
% Computational Cost: add. (301->143), mult. (310->192), div. (0->0), fcn. (168->8), ass. (0->70)
t79 = qJ(3) + qJ(4);
t74 = cos(t79);
t98 = t74 / 0.2e1;
t80 = qJ(1) + qJ(2);
t75 = cos(t80);
t97 = -t75 / 0.2e1;
t83 = cos(qJ(3));
t96 = t83 / 0.2e1;
t84 = cos(qJ(1));
t95 = t84 / 0.2e1;
t71 = V_base(6) + qJD(1);
t94 = pkin(2) * t71;
t81 = sin(qJ(3));
t93 = pkin(3) * t81;
t92 = pkin(3) * t83;
t82 = sin(qJ(1));
t91 = Icges(2,4) * t82;
t73 = sin(t80);
t90 = Icges(3,4) * t73;
t89 = Icges(3,4) * t75;
t88 = Icges(4,4) * t81;
t72 = sin(t79);
t87 = Icges(5,4) * t72;
t86 = V_base(6) * pkin(1) + V_base(2);
t70 = V_base(6) + qJD(3);
t77 = Icges(2,4) * t84;
t76 = Icges(4,4) * t83;
t69 = qJD(2) + t71;
t68 = qJD(4) + t70;
t67 = Icges(5,4) * t74;
t66 = t84 * rSges(2,1) - t82 * rSges(2,2);
t65 = t83 * rSges(4,1) - t81 * rSges(4,2);
t64 = t82 * rSges(2,1) + t84 * rSges(2,2);
t63 = t81 * rSges(4,1) + t83 * rSges(4,2);
t62 = Icges(2,1) * t84 - t91;
t61 = Icges(2,1) * t82 + t77;
t60 = Icges(4,1) * t83 - t88;
t59 = Icges(4,1) * t81 + t76;
t58 = -Icges(2,2) * t82 + t77;
t57 = Icges(2,2) * t84 + t91;
t56 = -Icges(4,2) * t81 + t76;
t55 = Icges(4,2) * t83 + t88;
t50 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t49 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t48 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t47 = -t75 * rSges(3,1) + t73 * rSges(3,2);
t46 = t74 * rSges(5,1) - t72 * rSges(5,2);
t45 = -t73 * rSges(3,1) - t75 * rSges(3,2);
t44 = t72 * rSges(5,1) + t74 * rSges(5,2);
t43 = -Icges(3,1) * t75 + t90;
t42 = -Icges(3,1) * t73 - t89;
t41 = Icges(5,1) * t74 - t87;
t40 = Icges(5,1) * t72 + t67;
t39 = Icges(3,2) * t73 - t89;
t38 = -Icges(3,2) * t75 - t90;
t37 = -Icges(5,2) * t72 + t67;
t36 = Icges(5,2) * t74 + t87;
t31 = V_base(5) * rSges(2,3) - t71 * t64 + V_base(1);
t30 = V_base(5) * rSges(4,3) - t70 * t63 + V_base(1);
t29 = -V_base(4) * rSges(2,3) + t71 * t66 + V_base(2);
t28 = -V_base(4) * rSges(4,3) + t70 * t65 + t86;
t27 = V_base(4) * t64 - V_base(5) * t66 + V_base(3);
t26 = V_base(4) * t63 + V_base(3) + (-pkin(1) - t65) * V_base(5);
t25 = V_base(5) * rSges(3,3) - t69 * t45 - t82 * t94 + V_base(1);
t24 = V_base(5) * rSges(5,3) - t68 * t44 - t70 * t93 + V_base(1);
t23 = -V_base(4) * rSges(3,3) + t69 * t47 + t84 * t94 + V_base(2);
t22 = -V_base(4) * rSges(5,3) + t68 * t46 + t70 * t92 + t86;
t21 = V_base(4) * t45 - V_base(5) * t47 + V_base(3) + (t82 * V_base(4) - t84 * V_base(5)) * pkin(2);
t20 = V_base(3) + (t44 + t93) * V_base(4) + (-pkin(1) - t46 - t92) * V_base(5);
t1 = m(1) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) / 0.2e1 + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + m(2) * (t27 ^ 2 + t29 ^ 2 + t31 ^ 2) / 0.2e1 + Icges(2,3) * t71 ^ 2 / 0.2e1 + m(3) * (t21 ^ 2 + t23 ^ 2 + t25 ^ 2) / 0.2e1 + Icges(3,3) * t69 ^ 2 / 0.2e1 + m(4) * (t26 ^ 2 + t28 ^ 2 + t30 ^ 2) / 0.2e1 + Icges(4,3) * t70 ^ 2 / 0.2e1 + m(5) * (t20 ^ 2 + t22 ^ 2 + t24 ^ 2) / 0.2e1 + Icges(5,3) * t68 ^ 2 / 0.2e1 + (Icges(1,6) * V_base(6) + (Icges(5,5) * t72 + Icges(5,6) * t74) * t68 + (-Icges(3,5) * t73 - Icges(3,6) * t75) * t69 + (Icges(4,5) * t81 + Icges(4,6) * t83) * t70 + (Icges(2,5) * t82 + Icges(2,6) * t84) * t71 + (Icges(1,2) / 0.2e1 + t57 * t95 + t82 * t61 / 0.2e1 + t38 * t97 - t73 * t42 / 0.2e1 + t55 * t96 + t81 * t59 / 0.2e1 + t36 * t98 + t72 * t40 / 0.2e1) * V_base(5)) * V_base(5) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + (Icges(2,5) * t84 - Icges(2,6) * t82) * t71 + (-Icges(3,5) * t75 + Icges(3,6) * t73) * t69 + (Icges(4,5) * t83 - Icges(4,6) * t81) * t70 + (Icges(5,5) * t74 - Icges(5,6) * t72) * t68 + ((t58 + t61) * t84 + (t56 + t59) * t83 + (-t57 + t62) * t82 + (-t55 + t60) * t81 + (-t39 - t42) * t75 + (t37 + t40) * t74 + (t38 - t43) * t73 + (-t36 + t41) * t72) * V_base(5) / 0.2e1 + (Icges(1,1) / 0.2e1 - t82 * t58 / 0.2e1 + t62 * t95 + t73 * t39 / 0.2e1 + t43 * t97 - t81 * t56 / 0.2e1 + t60 * t96 - t72 * t37 / 0.2e1 + t41 * t98) * V_base(4)) * V_base(4);
T = t1;
