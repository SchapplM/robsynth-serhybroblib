% Calculate joint inertia matrix with Newton Euler for
% fourbarprisIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
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
% Datum: 2020-05-07 09:59
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = fourbarprisIC_inertiaJ_snew_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisIC_inertiaJ_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisIC_inertiaJ_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisIC_inertiaJ_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisIC_inertiaJ_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisIC_inertiaJ_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:59:33
% EndTime: 2020-05-07 09:59:33
% DurationCPUTime: 0.03s
% Computational Cost: add. (18->10), mult. (22->14), div. (8->4), fcn. (22->4), ass. (0->7)
t22 = sin(qJ(3));
t23 = sin(qJ(1));
t24 = cos(qJ(3));
t25 = cos(qJ(1));
t26 = (-t24 * t22 - t25 * t23) / (pkin(3) + qJ(2)) / (-t24 ^ 2 + t25 ^ 2);
t18 = t22 * t25 - t23 * t24;
t1 = [m(3) + 0.1e1 / pkin(2) ^ 2 / t18 ^ 2 * Ifges(4,3) + ((2 * mrSges(3,1)) + (Ifges(3,2) + Ifges(2,3) + (m(3) * qJ(2) + (2 * mrSges(3,3))) * qJ(2)) * t26) * t26;];
Mq = t1;
