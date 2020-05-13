% Calculate joint inertia matrix for
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
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
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbarprisDE1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_inertiaJ_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_inertiaJ_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE1_inertiaJ_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisDE1_inertiaJ_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisDE1_inertiaJ_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:09:49
% EndTime: 2020-05-07 09:09:49
% DurationCPUTime: 0.15s
% Computational Cost: add. (105->18), mult. (58->17), div. (25->7), fcn. (2->2), ass. (0->8)
t15 = (qJ(1) + pkin(3));
t18 = -pkin(2) + t15;
t19 = -pkin(2) - t15;
t25 = (0.1e1 / (pkin(1) + t19) / (pkin(1) - t18) / (pkin(1) - t19) / (pkin(1) + t18));
t17 = (t15 ^ 2);
t16 = (qJ(1) ^ 2);
t2 = pkin(1) ^ 2 - pkin(2) ^ 2 - t16 - (2 * qJ(1) + pkin(3)) * pkin(3);
t1 = [2 * t2 / t15 * (-1 / t25) ^ (-0.1e1 / 0.2e1) * mrSges(3,1) + m(3) + (-4 * Ifges(4,3) * t17 + (-m(3) * t16 - 2 * mrSges(3,3) * qJ(1) - Ifges(3,2) - Ifges(2,3)) * t2 ^ 2 / t17) * t25;];
%% Postprocessing: Reshape Output
% From vec2symmat_1_matlab.m
res = [t1(1);];
Mq = res;
