% Calculate joint inertia matrix for
% fourbar2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
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
% Datum: 2020-04-24 20:37
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar2IC_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(2,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2IC_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2IC_inertiaJ_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2IC_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar2IC_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar2IC_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:37:51
% EndTime: 2020-04-24 20:37:51
% DurationCPUTime: 0.04s
% Computational Cost: add. (44->19), mult. (51->28), div. (7->3), fcn. (17->4), ass. (0->11)
t29 = pkin(2) * m(3);
t28 = rSges(3,1) * cos(qJ(2));
t24 = -qJ(3) + qJ(1);
t17 = sin(qJ(2) + t24);
t27 = (-pkin(1) * t17 + pkin(2) * sin(t24)) / t17;
t26 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t25 = t26 * m(3) + Icges(3,3);
t23 = 0.1e1 / pkin(1);
t19 = sin(qJ(2));
t16 = rSges(3,2) * t19 * t29;
t1 = [0.2e1 * t16 + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3) + t19 ^ 2 / t17 ^ 2 * ((rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3)) + ((pkin(2) - 0.2e1 * t28) * pkin(2) + t26) * m(3) + (0.2e1 * (-t28 * t29 + t16 + t25) * t23 + t23 ^ 2 * t25 * t27) * t27;];
Mq = t1;
