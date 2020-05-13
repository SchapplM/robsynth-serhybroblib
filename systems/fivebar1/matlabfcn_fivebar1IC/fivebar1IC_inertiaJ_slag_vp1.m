% Calculate joint inertia matrix for
% fivebar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fivebar1IC_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1IC_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1IC_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fivebar1IC_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:10
% EndTime: 2020-04-27 06:19:10
% DurationCPUTime: 0.11s
% Computational Cost: add. (212->49), mult. (312->72), div. (34->7), fcn. (96->11), ass. (0->36)
t52 = sin(qJ(2));
t81 = pkin(2) * t52;
t51 = sin(qJ(4));
t80 = pkin(3) * t51;
t71 = qJ(3) + qJ(4);
t72 = qJ(1) + qJ(2);
t75 = cos(0.2e1 * t71) - cos(0.2e1 * t72);
t39 = 0.1e1 / t75;
t55 = 0.2e1 * qJ(3);
t56 = 0.2e1 * qJ(1);
t79 = (t75 * pkin(4) + (cos(qJ(4) + t55) - cos(qJ(4) + t56 + 0.2e1 * qJ(2))) * pkin(3)) * t39;
t78 = (t75 * pkin(5) + (-cos(qJ(2) + 0.2e1 * qJ(4) + t55) + cos(qJ(2) + t56)) * pkin(2)) * t39;
t74 = rSges(5,1) ^ 2 + rSges(5,2) ^ 2;
t40 = t74 * m(5) + Icges(5,3);
t77 = t40 / pkin(4) ^ 2;
t73 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t41 = t73 * m(3) + Icges(3,3);
t76 = t41 / pkin(5) ^ 2;
t70 = pkin(2) * rSges(3,1) * cos(qJ(2));
t69 = rSges(5,2) * t80;
t63 = 0.1e1 / pkin(4);
t68 = t63 * t79;
t61 = 0.1e1 / pkin(5);
t67 = t61 * t78;
t48 = m(3) * rSges(3,2) * t81;
t38 = -m(3) * t70 + t41 + t48;
t47 = pkin(3) * m(5) * rSges(5,1) * cos(qJ(4));
t37 = -m(5) * t69 + t40 + t47;
t66 = pkin(2) ^ 2;
t65 = pkin(3) ^ 2;
t46 = -sin(t71 - t72);
t43 = 0.1e1 / t46 ^ 2;
t42 = 0.1e1 / t46;
t34 = -t41 * t67 + t38;
t33 = -t40 * t68 + t37;
t1 = [0.2e1 * t48 + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3) + t66 * t52 ^ 2 * t43 * t77 - (t34 + t38) * t67 + (t66 - 0.2e1 * t70 + t73) * m(3), (t34 * t61 * t80 + (t37 * t63 - t77 * t79) * t81) * t42; (t33 * t63 * t81 + (t38 * t61 - t76 * t78) * t80) * t42, t65 * t51 ^ 2 * t43 * t76 + 0.2e1 * t47 + (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3) + Icges(5,3) - (t33 + t37) * t68 + (t65 - 0.2e1 * t69 + t74) * m(5);];
Mq = t1;
