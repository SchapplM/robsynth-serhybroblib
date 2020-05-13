% Calculate joint inertia matrix for
% fourbar2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
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
% Datum: 2020-04-24 20:17
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar2TE_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(2,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2TE_inertiaJ_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2TE_inertiaJ_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2TE_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar2TE_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar2TE_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:17:25
% EndTime: 2020-04-24 20:17:25
% DurationCPUTime: 0.03s
% Computational Cost: add. (6->6), mult. (8->8), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + m(3) * pkin(2) ^ 2 + Icges(4,3) + Icges(2,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t1(1);];
Mq = res;
