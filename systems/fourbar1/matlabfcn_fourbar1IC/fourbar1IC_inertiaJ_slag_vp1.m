% Calculate joint inertia matrix for
% fourbar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
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
% Mq [1x1]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:15
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1IC_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1IC_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1IC_inertiaJ_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1IC_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1IC_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar1IC_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:15:50
% EndTime: 2020-04-24 20:15:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (44->19), mult. (54->30), div. (8->4), fcn. (17->4), ass. (0->11)
t27 = -qJ(3) + qJ(1);
t18 = sin(qJ(2) + t27);
t30 = (-pkin(3) * t18 + pkin(2) * sin(t27)) / t18;
t29 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t28 = t29 * m(3) + Icges(3,3);
t26 = rSges(3,1) * pkin(2) * cos(qJ(2));
t25 = pkin(2) ^ 2;
t24 = 0.1e1 / pkin(3);
t20 = sin(qJ(2));
t17 = m(3) * pkin(2) * rSges(3,2) * t20;
t1 = [0.2e1 * t17 + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3) + t20 ^ 2 * t25 / pkin(4) ^ 2 / t18 ^ 2 * ((rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3)) + (t25 - 0.2e1 * t26 + t29) * m(3) + (0.2e1 * (-m(3) * t26 + t17 + t28) * t24 + t24 ^ 2 * t28 * t30) * t30;];
Mq = t1;
