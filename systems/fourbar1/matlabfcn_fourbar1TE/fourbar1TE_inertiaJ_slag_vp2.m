% Calculate joint inertia matrix for
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [1x1]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbar1TE_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_inertiaJ_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_inertiaJ_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1TE_inertiaJ_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbar1TE_inertiaJ_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbar1TE_inertiaJ_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:49:12
% EndTime: 2020-04-24 19:49:12
% DurationCPUTime: 0.08s
% Computational Cost: add. (220->24), mult. (298->32), div. (14->5), fcn. (74->4), ass. (0->19)
t36 = pkin(3) ^ 2;
t31 = pkin(1) * cos(qJ(1));
t30 = pkin(1) ^ 2 - 0.2e1 * pkin(2) * t31;
t33 = -pkin(3) + pkin(4);
t34 = -pkin(3) - pkin(4);
t35 = ((pkin(2) - t34) * (pkin(2) + t34) + t30) * ((pkin(2) - t33) * (pkin(2) + t33) + t30);
t32 = pkin(1) * sin(qJ(1));
t29 = -pkin(4) ^ 2 + t36;
t28 = pkin(2) * t32;
t14 = pkin(2) ^ 2 + t30;
t25 = sqrt(-t35);
t15 = -pkin(2) + t31;
t12 = t14 + t29;
t5 = pkin(2) * t15 * t25;
t4 = t5 + (t14 - t29) * t28;
t3 = -t12 * t28 + t5;
t2 = (t12 * t15 + t25 * t32) * pkin(2);
t1 = t3 ^ 2;
t6 = [Ifges(2,3) + (m(3) * (t1 / 0.2e1 + t2 ^ 2 / 0.2e1) / t36 / 0.2e1 + (-t4 ^ 2 * Ifges(3,3) - t1 * Ifges(4,3)) / t35 + (-mrSges(3,1) * t2 + mrSges(3,2) * t3) / pkin(3) * t4 / t25) / t14 ^ 2;];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t6(1);];
Mq = res;
